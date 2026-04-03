/**
 * dipper_api.hpp
 *
 * Public C++ library API for DIPPER — phylogenetic placement of aligned MSA
 * sequences onto a backbone tree using GPU-accelerated k-closest placement.
 *
 * This header is the single include required to embed DIPPER as a library.
 * Link your target against the same object files that build the `dipper`
 * executable (see CMakeLists.txt) and add src/ to your include path.
 *
 * Two workflows are supported:
 *
 *   WITH backbone tree                 WITHOUT backbone tree
 *   ─────────────────────────────────  ─────────────────────────────────────
 *   1. loadBackboneTree()              (skip — no backbone needed)
 *   2. loadSequences()                 1. loadSequences()
 *   3. initializeGPU()                 2. initializeGPU()
 *   4. placeNextSequence() × N         3. placeNextSequence() × (N-2)
 *   5. writeTree()                     4. writeTree()
 *
 *   When no backbone tree is given, the first two sequences in the input seed
 *   an initial two-node tree; all remaining sequences are placed one at a time.
 *
 * Example
 * -------
 *   See examples/msa_placement_example.cpp.
 */

#pragma once

#include "fourBitCompressor.hpp"
#include "mash_placement.cuh"
#include "tree.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * PlacementResult
 *
 * Returned by DipperSession::placeNextSequence().
 *
 * Describes the edge split that occurred when one query sequence was placed
 * onto the tree.  Before placement, edge (splitNodeA ↔ splitNodeB) existed.
 * After placement, that edge is replaced by three new edges:
 *
 *   splitNodeA ──[lenA]── internalNodeIdx ──[lenB]── splitNodeB
 *                                 └──[lenTip]── tipNodeIdx  (query)
 *
 * The node indices are DIPPER's internal numbering.  Map them to sequence
 * names via DipperSession::sortedNames()[idx].
 *
 * For integration with an external node-pointer tree (N nodes before this
 * placement), create:
 *   • external node N+1 corresponding to internalNodeIdx
 *   • external node N+2 corresponding to tipNodeIdx
 * and connect them according to (splitNodeA, splitNodeB, lenA, lenB, lenTip).
 */
struct PlacementResult {
    bool more = false;      ///< true if more queries remain to be placed

    // ---- Split edge (DIPPER internal indices) --------------------------------
    int    splitNodeA = -1; ///< DIPPER idx: "from" end of the edge that was split
    int    splitNodeB = -1; ///< DIPPER idx: "to"   end of the edge that was split

    // ---- Split edge (sequential labels) -------------------------------------
    int    splitLabelA = -1; ///< sequential label of splitNodeA (use session.nodeLabel())
    int    splitLabelB = -1; ///< sequential label of splitNodeB

    // ---- New nodes (DIPPER internal indices) --------------------------------
    int    internalNodeIdx = -1; ///< DIPPER idx of the new internal node
    int    tipNodeIdx      = -1; ///< DIPPER idx of the new tip (query) node

    // ---- New nodes (sequential labels) --------------------------------------
    /// Sequential label of the new internal node.
    /// If N nodes existed before this placement, internalLabel = N.
    int    internalLabel = -1;
    /// Sequential label of the new tip (query) node.
    /// If N nodes existed before this placement, tipLabel = N+1.
    int    tipLabel = -1;

    std::string tipName;    ///< sequence name of the placed query

    // ---- Branch lengths (all >= 0) ------------------------------------------
    double lenA   = 0.0;    ///< splitNodeA  → internalNodeIdx
    double lenB   = 0.0;    ///< splitNodeB  → internalNodeIdx
    double lenTip = 0.0;    ///< internalNodeIdx → tipNodeIdx (query)
};

/**
 * DipperSession
 *
 * Encapsulates one complete DIPPER MSA-placement session.  Each session owns
 * its GPU device arrays independently of the static singletons used by the
 * CLI, so multiple sessions can be constructed in the same process (serially).
 *
 * The class is non-copyable because it owns GPU resources.
 */
class DipperSession {
public:
    DipperSession()  = default;
    ~DipperSession() { _cleanup(); }

    DipperSession(const DipperSession&)            = delete;
    DipperSession& operator=(const DipperSession&) = delete;

    // =========================================================================
    // Step 1 (optional) — Backbone tree
    // =========================================================================

    /**
     * Parse @p newickStr as the backbone tree in Newick format.
     *
     * This step is OPTIONAL.  If omitted, the first two sequences supplied to
     * loadSequences() will seed the initial tree, and all remaining sequences
     * will be treated as queries.
     *
     * @param newickStr   A complete Newick string (single line, terminated by ';').
     *
     * Must be called before loadSequences() when a backbone tree is desired.
     */
    void loadBackboneTree(const std::string& newickStr)
    {
        _requireOneOf({State::Empty}, "loadBackboneTree");
        // Store the raw Newick string.  The Tree object is constructed (with
        // the correct totalLeaves count) later in loadSequences(), once the
        // full sequence set is known.  A temporary parse here is only used to
        // read m_backboneSize for the diagnostic message.
        m_backboneNewick = newickStr;
        delete m_tree;
        m_tree         = new Tree(newickStr);   // temporary; rebuilt in loadSequences
        m_backboneSize = m_tree->m_numLeaves;
        m_state        = State::TreeLoaded;
        std::cerr << "[DIPPER] Backbone tree stored: " << m_backboneSize
                  << " leaves, root = " << m_tree->root->name << "\n";
    }

    /**
     * Read the first line of @p filePath and parse it as a Newick backbone tree.
     *
     * @return true on success, false if the file cannot be opened.
     */
    bool loadBackboneTreeFromFile(const std::string& filePath)
    {
        _requireOneOf({State::Empty}, "loadBackboneTreeFromFile");
        std::ifstream f(filePath);
        if (!f) {
            std::cerr << "[DIPPER] ERROR: cannot open backbone tree file: " << filePath << "\n";
            return false;
        }
        std::string newick;
        std::getline(f, newick);
        loadBackboneTree(newick);
        return true;
    }

    // =========================================================================
    // Step 2 — Aligned sequences
    // =========================================================================

    /**
     * Supply all aligned sequences as pre-loaded in-memory vectors.
     *
     * WITH backbone tree
     * ------------------
     * Each sequence is matched against the backbone tree by name:
     *   • If the name exists in the backbone, the sequence fills its tree-index
     *     slot (used only for distance computation).
     *   • If the name is absent, the sequence is appended as a query and will
     *     be placed during the placeNextSequence() calls.
     *
     * WITHOUT backbone tree
     * ---------------------
     * All sequences are accepted as-is in the order provided.  The first two
     * sequences seed the initial two-node tree; all remaining sequences are
     * placed as queries.  At least 3 sequences are required.
     *
     * @param sequences  Aligned sequence strings (all the same length for MSA).
     * @param names      Sequence names corresponding to each entry in @p sequences.
     */
    void loadSequences(const std::vector<std::string>& sequences,
                       const std::vector<std::string>& names)
    {
        _requireOneOf({State::Empty, State::TreeLoaded}, "loadSequences");
        if (sequences.size() != names.size())
            throw std::invalid_argument(
                "[DIPPER] loadSequences: sequences and names must have the same size.");

        const bool hasBackbone = (m_tree != nullptr);
        if (hasBackbone && sequences.empty())
            throw std::invalid_argument(
                "[DIPPER] loadSequences: at least one sequence is required.");
        if (!hasBackbone && sequences.size() < 3)
            throw std::invalid_argument(
                "[DIPPER] loadSequences: at least 3 sequences are required "
                "when no backbone tree is provided (first two seed the initial tree).");

        m_rawSeqs      = sequences;
        m_rawNames     = names;
        m_numSequences = sequences.size();

        if (hasBackbone) {
            // Rebuild the tree with totalLeaves = m_numSequences so that
            // internal node indices are offset above (m_numSequences - 1) and
            // cannot collide with query sequence slots.  This matches the CLI
            // behaviour (tree_generation.cu: new Tree(newick, namesDump.size())).
            delete m_tree;
            m_tree = new Tree(m_backboneNewick, m_numSequences);
            m_backboneSize = m_tree->m_numLeaves;
            _buildIdMapWithBackbone();
            _assignPreorderLabels();   // label all backbone nodes 0..K-1
            std::cerr << "[DIPPER] Sequences loaded: " << m_numSequences
                      << " total (" << m_backboneSize << " backbone, "
                      << (m_numSequences - m_backboneSize) << " queries).\n" << std::flush;
        } else {
            _buildIdMapNoBackbone();  // seeds get labels 0 and 1
            std::cerr << "[DIPPER] Sequences loaded: " << m_numSequences
                      << " total (no backbone — first 2 will seed the tree, "
                      << (m_numSequences - 2) << " queries)." << std::flush;
        }
        m_state = State::SeqsLoaded;
    }

    // =========================================================================
    // Step 3 — GPU initialization
    // =========================================================================

    /**
     * Compress all sequences with 4-bit encoding, upload them to GPU memory,
     * and prepare the initial tree topology on the GPU.
     *
     * With a backbone tree: loads the provided tree topology onto the GPU.
     * Without a backbone tree: computes the distance between the first two
     *   sequences and builds the initial two-node tree on the GPU.
     *
     * The @p params object is updated in-place:
     *   • params.in  is forced to "m" (MSA / aligned input).
     *   • params.out is forced to "t" (tree output).
     *   • params.totalNumSeqs and params.backboneSize are set automatically.
     *
     * All other fields (range, isProtein, threshold, distanceType, etc.) are
     * respected as supplied by the caller.
     *
     * Must be called after loadSequences().
     */
    void initializeGPU(MashPlacement::Param& params)
    {
        _requireOneOf({State::SeqsLoaded}, "initializeGPU");
        m_params            = &params;
        params.in           = "m";
        params.out          = "t";
        params.totalNumSeqs = m_numSequences;
        params.backboneSize = m_backboneSize;   // 0 when no backbone

        std::cerr << "[DIPPER] Compressing sequences (4-bit)...\n";
        _compressSequences(params);

        std::cerr << "[DIPPER] Uploading sequences to GPU...\n";
        m_msaArrays.allocateDeviceArrays(
            m_fourBitSeqs, m_seqLengths, m_numSequences, params);

        const bool hasBackbone = (m_tree != nullptr);
        if (hasBackbone) {
            std::cerr << "[DIPPER] Loading backbone tree topology onto GPU...\n";
            m_kpArrays.allocateDeviceArrays(m_numSequences, (int)m_backboneSize);
            m_kpArrays.initializeDeviceArrays(m_tree);
            m_kpArrays.beginIterativePlacement();
            std::cerr << "[DIPPER] GPU initialized. Ready to place "
                      << (m_numSequences - m_backboneSize) << " query sequence(s).\n";
        } else {
            std::cerr << "[DIPPER] Building initial two-sequence tree on GPU...\n";
            // backboneSize == -1 signals scratch mode to allocateDeviceArrays.
            m_kpArrays.allocateDeviceArrays(m_numSequences, /*backboneSize=*/-1);
            m_kpArrays.beginIterativePlacementFromScratch(
                params, m_mashArrays, m_matrixReader, m_msaArrays);
            std::cerr << "[DIPPER] GPU initialized. Ready to place "
                      << (m_numSequences - 2) << " query sequence(s).\n";
        }

        m_state = State::GPUReady;
    }

    // =========================================================================
    // Step 4 — Place one query sequence at a time
    // =========================================================================

    /**
     * Place the next query sequence onto the current tree and return the
     * edge-split details needed to update an external node-pointer tree.
     *
     * Each call places exactly one sequence.  Calls must be made in order
     * after initializeGPU().
     *
     * With backbone:    queries start from index backboneSize.
     * Without backbone: queries start from index 2 (first two seeded the tree).
     *
     * @return A PlacementResult whose fields describe:
     *   • more            — whether additional queries remain.
     *   • splitNodeA/B    — DIPPER indices of the edge endpoints that were split.
     *   • internalNodeIdx — DIPPER index of the newly inserted internal node.
     *   • tipNodeIdx      — DIPPER index of the newly inserted tip (query).
     *   • tipName         — sequence name of the placed query.
     *   • lenA/lenB/lenTip — branch lengths for the three new edges.
     *
     * @throws std::logic_error if all queries have already been placed.
     */
    PlacementResult placeNextSequence()
    {
        _requireOneOf({State::GPUReady}, "placeNextSequence");
        if (m_placementPending)
            throw std::logic_error(
                "[DIPPER] placeNextSequence: a two-phase placement is in progress. "
                "Call commitPlacement() first.");
        if (m_kpArrays.m_nextQueryIdx >= (int)m_numSequences)
            throw std::logic_error(
                "[DIPPER] placeNextSequence: all query sequences have already been placed.");

        const int seqIdx    = m_kpArrays.m_nextQueryIdx;
        const int seedCount = (m_tree != nullptr) ? (int)m_backboneSize : 2;
        const int queryNum  = seqIdx - seedCount + 1;
        const int totalQ    = (int)m_numSequences - seedCount;

        std::cerr << "[DIPPER] Placing query " << queryNum << " / " << totalQ
                  << " (global idx " << seqIdx
                  << ", name: " << m_sortedNames[seqIdx] << ")...\n";

        MashPlacement::KPlacementDeviceArrays::PlacementInfo info{};
        bool more = m_kpArrays.placeSingleQuery(
            seqIdx, *m_params, m_mashArrays, m_matrixReader, m_msaArrays, &info);

        // Assign sequential labels to the two new nodes and register them.
        const int iLabel = m_nextLabel;       // new internal node label (N)
        const int tLabel = m_nextLabel + 1;   // new tip node label      (N+1)
        m_dipperToLabel[info.internalNodeIdx] = iLabel;
        m_dipperToLabel[info.tipNodeIdx]      = tLabel;
        m_labelToDipper[iLabel]               = info.internalNodeIdx;
        m_labelToDipper[tLabel]               = info.tipNodeIdx;
        m_nextLabel += 2;

        PlacementResult result;
        result.more            = more;
        result.splitNodeA      = info.splitNodeA;
        result.splitNodeB      = info.splitNodeB;
        result.splitLabelA     = nodeLabel(info.splitNodeA);
        result.splitLabelB     = nodeLabel(info.splitNodeB);
        result.internalNodeIdx = info.internalNodeIdx;
        result.tipNodeIdx      = info.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_sortedNames[seqIdx];
        result.lenA            = info.lenA;
        result.lenB            = info.lenB;
        result.lenTip          = info.lenTip;

        if (!more) {
            m_state = State::Done;
            std::cerr << "[DIPPER] All query sequences placed.\n";
        }
        return result;
    }

    // =========================================================================
    // Step 4 — Two-phase placement: findNextPlacement / commitPlacement
    // =========================================================================

    /**
     * Phase 1 of two-phase placement.
     *
     * Computes DIPPER's optimal placement for the next query sequence and
     * returns the edge-split details, but does NOT insert the tip into the
     * GPU tree yet.  Call commitPlacement() to finalise.
     *
     * The returned PlacementResult reflects the edge that DIPPER recommends.
     * The caller may override this choice in commitPlacement() by supplying
     * different node labels.
     *
     * @throws std::logic_error if commitPlacement() has not been called after
     *         a previous findNextPlacement(), or if all queries are placed.
     */
    PlacementResult findNextPlacement()
    {
        _requireOneOf({State::GPUReady}, "findNextPlacement");
        if (m_placementPending)
            throw std::logic_error(
                "[DIPPER] findNextPlacement: a placement is already pending. "
                "Call commitPlacement() first.");
        if (m_kpArrays.m_nextQueryIdx >= (int)m_numSequences)
            throw std::logic_error(
                "[DIPPER] findNextPlacement: all query sequences have already been placed.");

        const int seqIdx    = m_kpArrays.m_nextQueryIdx;
        const int seedCount = (m_tree != nullptr) ? (int)m_backboneSize : 2;
        const int queryNum  = seqIdx - seedCount + 1;
        const int totalQ    = (int)m_numSequences - seedCount;

        std::cerr << "[DIPPER] Finding placement for query " << queryNum << " / " << totalQ
                  << " (global idx " << seqIdx
                  << ", name: " << m_sortedNames[seqIdx] << ")...\n";

        MashPlacement::KPlacementDeviceArrays::PlacementInfo info =
            m_kpArrays.findBestEdge(seqIdx, *m_params, m_mashArrays, m_matrixReader, m_msaArrays);

        // Pre-assign sequential labels (fixed regardless of which edge is chosen).
        const int iLabel = m_nextLabel;
        const int tLabel = m_nextLabel + 1;

        PlacementResult result;
        result.more            = (seqIdx + 1 < (int)m_numSequences);
        result.splitNodeA      = info.splitNodeA;
        result.splitNodeB      = info.splitNodeB;
        result.splitLabelA     = nodeLabel(info.splitNodeA);
        result.splitLabelB     = nodeLabel(info.splitNodeB);
        result.internalNodeIdx = info.internalNodeIdx;
        result.tipNodeIdx      = info.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_sortedNames[seqIdx];
        result.lenA            = info.lenA;
        result.lenB            = info.lenB;
        result.lenTip          = info.lenTip;

        m_pendingResult    = result;
        m_placementPending = true;
        return result;
    }

    /**
     * Phase 2 of two-phase placement.
     *
     * Inserts the pending query sequence into the GPU tree on the edge whose
     * two endpoints carry the given sequential labels.  The labels may be the
     * same pair returned by findNextPlacement() (confirming DIPPER's choice)
     * or a different pair (overriding it — the chosen edge must exist in the
     * current tree).
     *
     * @param labelA  Sequential label of one endpoint of the target edge.
     * @param labelB  Sequential label of the other endpoint.
     * @return        Final PlacementResult for this query, reflecting the
     *                actual edge used (branch lengths/split info may differ
     *                from findNextPlacement() if an override was requested).
     *
     * @throws std::logic_error if findNextPlacement() has not been called, or
     *         if the edge specified by (labelA, labelB) cannot be found.
     */
    PlacementResult commitPlacement(int labelA, int labelB)
    {
        _requireOneOf({State::GPUReady}, "commitPlacement");
        if (!m_placementPending)
            throw std::logic_error(
                "[DIPPER] commitPlacement: no pending placement. "
                "Call findNextPlacement() first.");

        const int seqIdx = m_kpArrays.m_pendingSeqIdx;

        // Map sequential labels → DIPPER node indices.
        auto itA = m_labelToDipper.find(labelA);
        auto itB = m_labelToDipper.find(labelB);
        if (itA == m_labelToDipper.end() || itB == m_labelToDipper.end())
            throw std::logic_error(
                "[DIPPER] commitPlacement: one or both labels are unknown: "
                + std::to_string(labelA) + ", " + std::to_string(labelB));

        const int dipperA = itA->second;
        const int dipperB = itB->second;

        // Compare against the split edge already stored in m_pendingResult
        // (populated by findBestEdge → getEdgeInfo, no extra GPU reads needed).
        const int defaultDipperA = m_pendingResult.splitNodeA;
        const int defaultDipperB = m_pendingResult.splitNodeB;

        const bool isOverride =
            !((dipperA == defaultDipperA && dipperB == defaultDipperB) ||
              (dipperA == defaultDipperB && dipperB == defaultDipperA));

        int overrideEid = -1;
        MashPlacement::KPlacementDeviceArrays::PlacementInfo finalInfo = {};

        if (isOverride) {
            overrideEid = m_kpArrays.findEdgeBetween(dipperA, dipperB);
            if (overrideEid < 0)
                throw std::logic_error(
                    "[DIPPER] commitPlacement: no edge found between DIPPER nodes "
                    + std::to_string(dipperA) + " and " + std::to_string(dipperB)
                    + " (labels " + std::to_string(labelA) + " / " + std::to_string(labelB) + ").");
            finalInfo = m_kpArrays.getEdgeInfo(overrideEid, seqIdx);
        } else {
            finalInfo.splitNodeA      = m_pendingResult.splitNodeA;
            finalInfo.splitNodeB      = m_pendingResult.splitNodeB;
            finalInfo.internalNodeIdx = m_pendingResult.internalNodeIdx;
            finalInfo.tipNodeIdx      = m_pendingResult.tipNodeIdx;
            finalInfo.lenA            = m_pendingResult.lenA;
            finalInfo.lenB            = m_pendingResult.lenB;
            finalInfo.lenTip          = m_pendingResult.lenTip;
        }

        // Commit the placement on the GPU.
        m_kpArrays.commitQuery(seqIdx, overrideEid);

        // Register the two new nodes in both label maps.
        const int iLabel = m_pendingResult.internalLabel;
        const int tLabel = m_pendingResult.tipLabel;
        m_dipperToLabel[finalInfo.internalNodeIdx] = iLabel;
        m_dipperToLabel[finalInfo.tipNodeIdx]      = tLabel;
        m_labelToDipper[iLabel]                    = finalInfo.internalNodeIdx;
        m_labelToDipper[tLabel]                    = finalInfo.tipNodeIdx;
        m_nextLabel += 2;

        const bool more = (m_kpArrays.m_nextQueryIdx < (int)m_numSequences);

        PlacementResult result;
        result.more            = more;
        result.splitNodeA      = finalInfo.splitNodeA;
        result.splitNodeB      = finalInfo.splitNodeB;
        result.splitLabelA     = nodeLabel(finalInfo.splitNodeA);
        result.splitLabelB     = nodeLabel(finalInfo.splitNodeB);
        result.internalNodeIdx = finalInfo.internalNodeIdx;
        result.tipNodeIdx      = finalInfo.tipNodeIdx;
        result.internalLabel   = iLabel;
        result.tipLabel        = tLabel;
        result.tipName         = m_pendingResult.tipName;
        result.lenA            = finalInfo.lenA;
        result.lenB            = finalInfo.lenB;
        result.lenTip          = finalInfo.lenTip;

        m_placementPending = false;
        m_pendingResult    = {};

        if (!more) {
            m_state = State::Done;
            std::cerr << "[DIPPER] All query sequences placed.\n";
        }
        return result;
    }

    // =========================================================================
    // Output
    // =========================================================================

    /**
     * Write the placement tree in Newick format to the file at @p outputPath.
     *
     * Can be called at any point after initializeGPU(), even before all
     * queries are placed, to inspect intermediate results.
     */
    void writeTree(const std::string& outputPath)
    {
        if (m_state < State::GPUReady)
            throw std::logic_error(
                "[DIPPER] writeTree: initializeGPU() must be called first.");
        std::ofstream ofs(outputPath);
        if (!ofs)
            throw std::runtime_error(
                "[DIPPER] writeTree: cannot open output file: " + outputPath);
        m_kpArrays.printTree(m_sortedNames, ofs);
        std::cerr << "[DIPPER] Tree written to " << outputPath << "\n";
    }

    /** Write the placement tree to an already-open output stream. */
    void writeTree(std::ofstream& ofs)
    {
        if (m_state < State::GPUReady)
            throw std::logic_error(
                "[DIPPER] writeTree: initializeGPU() must be called first.");
        m_kpArrays.printTree(m_sortedNames, ofs);
    }

    // =========================================================================
    // Convenience: place ALL remaining queries in one call
    // =========================================================================

    /**
     * Place every remaining query sequence in insertion order.
     *
     * Equivalent to: while (placeNextSequence()) {}
     */
    void placeAllSequences()
    {
        while (m_kpArrays.m_nextQueryIdx < (int)m_numSequences)
            (void)placeNextSequence();
    }

    // =========================================================================
    // Accessors
    // =========================================================================

    /** Number of leaves in the backbone tree (0 when no backbone was given). */
    size_t backboneSize()   const { return m_backboneSize; }

    /** True when no backbone tree was provided (from-scratch mode). */
    bool   hasBackbone()    const { return m_tree != nullptr; }

    /** Total number of sequences (backbone + query, or all sequences). */
    size_t numSequences()   const { return m_numSequences; }

    /** Total number of nodes in the backbone tree (leaves + internal nodes).
     *
     *  This is the correct starting value for an external N counter:
     *    - Without backbone: returns 2 (the two seed nodes).
     *    - With backbone:    returns allNodes.size() from the parsed Newick.
     *
     *  Each call to placeNextSequence() adds 2 nodes, so after k placements
     *  the external tree has backboneNodeCount() + 2*k nodes. */
    size_t backboneNodeCount() const
    {
        if (m_tree) return m_tree->allNodes.size();
        return 2;   // no-backbone mode: two seed nodes
    }

    /** Global index of the next sequence that placeNextSequence() will place. */
    int    nextQueryIndex() const { return m_kpArrays.m_nextQueryIdx; }

    /**
     * Convert a DIPPER internal node index to its sequential pre-order label.
     *
     * Valid for any DIPPER index that has appeared in PlacementResult
     * (splitNodeA, splitNodeB, internalNodeIdx, tipNodeIdx) or in the
     * backbone tree nodes.
     *
     * @return  The 0-based sequential label, or -1 if the index is unknown.
     */
    int nodeLabel(int dipperIdx) const
    {
        auto it = m_dipperToLabel.find(dipperIdx);
        return (it != m_dipperToLabel.end()) ? it->second : -1;
    }

    /** Names in sorted GPU order.
     *
     *  Layout (with backbone):
     *    [0 .. backboneSize-1]          backbone leaf names
     *    [backboneSize .. numSeqs-1]    query names (not yet placed)
     *    [numSeqs .. ]                  internal backbone node names (e.g. "node_1")
     *
     *  Accessing splitNodeA/splitNodeB from PlacementResult is always safe
     *  because _buildIdMapWithBackbone() ensures internal nodes are included. */
    const std::vector<std::string>& sortedNames() const { return m_sortedNames; }

private:
    // ---- Session state machine -----------------------------------------------
    enum class State { Empty, TreeLoaded, SeqsLoaded, GPUReady, Done };
    State m_state = State::Empty;

    void _requireOneOf(std::initializer_list<State> allowed, const char* caller) const
    {
        for (State s : allowed)
            if (m_state == s) return;
        throw std::logic_error(
            std::string("[DIPPER] ") + caller + "() called out of order "
            "(current state: " + std::to_string(static_cast<int>(m_state)) + ").");
    }

    // ---- Owned data ----------------------------------------------------------
    Tree*                        m_tree           = nullptr;  // null = no backbone
    std::string                  m_backboneNewick;            // stored for deferred rebuild
    std::vector<std::string>     m_rawSeqs;
    std::vector<std::string>     m_rawNames;
    std::vector<std::string>     m_sortedNames;
    std::unordered_map<int, int> m_idMap;
    size_t m_numSequences = 0;
    size_t m_backboneSize = 0;   // 0 when no backbone

    // ---- Sequential node labeling (pre-order DFS, 0-based) ------------------
    // Backbone nodes are labeled 0..K-1 at init; each placement adds two more.
    std::unordered_map<int, int> m_dipperToLabel;  // DIPPER idx → sequential label
    std::unordered_map<int, int> m_labelToDipper;  // sequential label → DIPPER idx (reverse)
    int m_nextLabel = 0;                           // next label to assign

    // ---- Two-phase placement state ------------------------------------------
    bool            m_placementPending = false;
    PlacementResult m_pendingResult    {};          // filled by findNextPlacement()

    uint64_t** m_fourBitSeqs = nullptr;
    uint64_t*  m_seqLengths  = nullptr;

    MashPlacement::Param* m_params = nullptr;

    // GPU device array objects (each session owns its own, not the CLI statics).
    MashPlacement::MSADeviceArrays        m_msaArrays    {};
    MashPlacement::KPlacementDeviceArrays m_kpArrays     {};
    // Unused for MSA input, but required by existing method signatures:
    MashPlacement::MashDeviceArrays       m_mashArrays   {};
    MashPlacement::MatrixReader           m_matrixReader {};

    // ---- Helpers -------------------------------------------------------------

    /**
     * Backbone path: backbone sequences fill tree-index slots; novel sequences
     * are appended as queries; internal backbone nodes are appended beyond the
     * query slots so that sortedNames[idx] is always valid for any DIPPER node
     * index that can appear in PlacementResult::splitNodeA / splitNodeB.
     *
     * Final layout of sortedNames:
     *   [0 .. backboneSize-1]   backbone leaf names (at their tree-assigned idx)
     *   [backboneSize .. numSeqs-1]  query names (in input order)
     *   [numSeqs .. ]           internal backbone node names (e.g. "node_1")
     */
    void _buildIdMapWithBackbone()
    {
        m_sortedNames.assign(m_backboneSize, "");
        m_idMap.clear();

        for (int i = 0; i < (int)m_numSequences; ++i) {
            auto it = m_tree->allNodes.find(m_rawNames[i]);
            if (it == m_tree->allNodes.end()) {
                // Novel query: append after backbone slots.
                m_sortedNames.push_back(m_rawNames[i]);
                m_idMap[i] = (int)m_sortedNames.size() - 1;
            } else {
                // Backbone sequence: fill its tree-index slot.
                m_sortedNames[it->second->idx] = m_rawNames[i];
                m_idMap[i] = it->second->idx;
            }
        }

        // Internal backbone nodes have idx >= m_numSequences (because the Tree
        // was constructed with totalLeaves = m_numSequences, pushing internal
        // node indices above the sequence slots).  Extend sortedNames so that
        // sortedNames[node->idx] gives the internal node's Newick name, making
        // PlacementResult::splitNodeA/B always safe to index.
        for (auto& [name, node] : m_tree->allNodes) {
            if (!node->children.empty()) {  // internal node
                const int nidx = node->idx;
                if (nidx >= (int)m_sortedNames.size())
                    m_sortedNames.resize(nidx + 1, "");
                m_sortedNames[nidx] = name;   // e.g. "node_1", "node_2", "node_3"
            }
        }
    }

    /**
     * Pre-order DFS of the backbone tree to assign sequential labels 0, 1, 2, …
     * to every node (root first, then each subtree in child order).
     * After this call, m_nextLabel == backboneNodeCount().
     */
    void _assignPreorderLabels()
    {
        m_dipperToLabel.clear();
        m_labelToDipper.clear();
        m_nextLabel = 0;
        std::function<void(Node*)> dfs = [&](Node* node) {
            if (!node) return;
            const int label = m_nextLabel++;
            m_dipperToLabel[node->idx] = label;
            m_labelToDipper[label]     = node->idx;
            for (Node* child : node->children)
                dfs(child);
        };
        dfs(m_tree->root);
    }

    /**
     * No-backbone path: sequences retain their input order (0, 1, 2, …).
     * Indices 0 and 1 seed the initial tree; the rest are queries.
     * Seed nodes get sequential labels 0 and 1.
     */
    void _buildIdMapNoBackbone()
    {
        m_backboneSize = 0;
        m_sortedNames  = m_rawNames;
        m_idMap.clear();
        for (int i = 0; i < (int)m_numSequences; ++i)
            m_idMap[i] = i;

        // Seed nodes 0 and 1 get the first two labels.
        m_dipperToLabel.clear();
        m_labelToDipper.clear();
        m_dipperToLabel[0] = 0;  m_labelToDipper[0] = 0;
        m_dipperToLabel[1] = 1;  m_labelToDipper[1] = 1;
        m_nextLabel = 2;
    }

    /** 4-bit compress all sequences in parallel using TBB. */
    void _compressSequences(const MashPlacement::Param& params)
    {
        m_fourBitSeqs = new uint64_t*[m_numSequences];
        m_seqLengths  = new uint64_t [m_numSequences];

        const bool rangeModify =
            (params.range.first > 0 || params.range.second > -1);

        tbb::parallel_for(
            tbb::blocked_range<int>(0, (int)m_numSequences),
            [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i < r.end(); ++i) {
                    int seqLen = (int)m_rawSeqs[i].size();
                    if (rangeModify) {
                        if (params.range.second > -1) seqLen = params.range.second + 1;
                        if (params.range.first  > 0)  seqLen -= params.range.first;
                    }
                    const uint64_t compSize = params.isProtein
                        ? (uint64_t)(seqLen + 7)  / 8
                        : (uint64_t)(seqLen + 15) / 16;
                    uint64_t* buf = new uint64_t[compSize];
                    fourBitCompressor(
                        m_rawSeqs[i], m_rawSeqs[i].size(), buf,
                        params.range.first, params.range.second, params.isProtein);
                    const int newId      = m_idMap[i];
                    m_seqLengths[newId]  = (uint64_t)seqLen;
                    m_fourBitSeqs[newId] = buf;
                }
            });
    }

    void _cleanup()
    {
        if (m_state >= State::GPUReady) {
            m_kpArrays.endIterativePlacement();
            m_msaArrays.deallocateDeviceArrays();
            m_kpArrays.deallocateDeviceArrays();
        }
        if (m_fourBitSeqs) {
            for (size_t i = 0; i < m_numSequences; ++i)
                delete[] m_fourBitSeqs[i];
            delete[] m_fourBitSeqs;
            m_fourBitSeqs = nullptr;
        }
        delete[] m_seqLengths;
        m_seqLengths = nullptr;
        delete m_tree;
        m_tree = nullptr;
    }
};
