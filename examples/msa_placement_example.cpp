/**
 * examples/msa_placement_example.cpp
 *
 * Demonstrates both usage modes of the DIPPER library API (DipperSession):
 *
 *   Mode A — WITH backbone tree (placement onto an existing tree)
 *     ./dipper_example_msa --backbone backbone.nwk sequences.fa output.nwk
 *
 *   Mode B — WITHOUT backbone tree (build tree from scratch, one seq at a time)
 *     ./dipper_example_msa sequences.fa output.nwk
 *
 * Input expectations
 * ------------------
 *   backbone.nwk  — single-line Newick file (Mode A only).
 *   sequences.fa  — FASTA file with aligned sequences (all the same length).
 *                   Mode A: must contain all backbone leaves plus query sequences.
 *                   Mode B: first two sequences seed the initial tree; all
 *                           remaining sequences are placed one at a time.
 *   output.nwk    — path where the resulting Newick tree is written.
 *
 * Placement API
 * -------------
 * Two approaches are available; choose whichever fits your integration:
 *
 *   One-phase (simple):
 *     PlacementResult r = session.placeNextSequence();
 *     // GPU tree updated immediately.  Use r.splitLabelA/B, r.lenA/B/Tip, etc.
 *
 *   Two-phase (for tight integration with an external node-pointer tree):
 *     PlacementResult proposal = session.findNextPlacement();
 *     // Inspect proposal: DIPPER suggests placing on the edge between
 *     //   proposal.splitLabelA  and  proposal.splitLabelB.
 *     // Here you can update your external tree and/or choose a different edge.
 *     PlacementResult final = session.commitPlacement(
 *         proposal.splitLabelA,   // confirm DIPPER's edge …
 *         proposal.splitLabelB);  // … or pass different labels to override.
 *     // GPU tree is updated only after commitPlacement().
 *
 * Build (from the examples/ directory):
 *   mkdir build && cd build
 *   cmake -DUSE_CUDA=ON ..
 *   make dipper_example_msa
 */

#include "dipper_api.hpp"   // DipperSession — the only DIPPER header needed

#include <zlib.h>
#include "kseq.h"           // lightweight kseq FASTA parser (must follow zlib.h)
KSEQ_INIT2(, gzFile, gzread)

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Minimal FASTA reader using kseq.
// ---------------------------------------------------------------------------
static bool readFasta(const std::string& path,
                      std::vector<std::string>& seqs,
                      std::vector<std::string>& names)
{
    gzFile f = gzopen(path.c_str(), "r");
    if (!f) {
        std::cerr << "ERROR: cannot open FASTA file: " << path << "\n";
        return false;
    }
    kseq_t* kseq = kseq_init(f);
    while (kseq_read(kseq) >= 0) {
        names.push_back(std::string(kseq->name.s, kseq->name.l));
        seqs.push_back(std::string(kseq->seq.s,  kseq->seq.l));
    }
    kseq_destroy(kseq);
    gzclose(f);
    return true;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // Parse arguments — accept both modes.
    std::string backbonePath, fastaPath, outputPath;
    bool withBackbone = false;

    if (argc == 4 && std::string(argv[1]) == "--backbone") {
        // Mode A: --backbone <backbone.nwk> <sequences.fa> <output.nwk>
        backbonePath = argv[2];
        fastaPath    = argv[3];
        outputPath   = argv[4];
        withBackbone = true;
    } else if (argc == 5 && std::string(argv[1]) == "--backbone") {
        backbonePath = argv[2];
        fastaPath    = argv[3];
        outputPath   = argv[4];
        withBackbone = true;
    } else if (argc == 3) {
        // Mode B: <sequences.fa> <output.nwk>
        fastaPath  = argv[1];
        outputPath = argv[2];
    } else {
        std::cerr << "Usage (with backbone):    " << argv[0]
                  << " --backbone <backbone.nwk> <sequences.fa> <output.nwk>\n"
                  << "Usage (without backbone): " << argv[0]
                  << " <sequences.fa> <output.nwk>\n";
        return 1;
    }

    // -------------------------------------------------------------------------
    // Configure DIPPER parameters.
    //
    // For MSA input, kmerSize / sketchSize / threshold are unused; set
    // distanceType=0 for Hamming distance on aligned sequences.
    // -------------------------------------------------------------------------
    MashPlacement::Param params(
        /*kmerSize=*/   21,
        /*sketchSize=*/ 1000,
        /*threshold=*/  0,
        /*distanceType*/0,
        /*in=*/         "m",   // forced to "m" by initializeGPU()
        /*out=*/        "t"
    );
    params.range     = {0, -1};   // use the full alignment length
    params.isProtein = false;     // DNA / RNA sequences
    params.distanceType = 2;

    // -------------------------------------------------------------------------
    // Read aligned sequences from FASTA.
    // -------------------------------------------------------------------------
    std::vector<std::string> seqs, names;
    if (!readFasta(fastaPath, seqs, names)) return 1;
    std::cerr << "Read " << seqs.size() << " sequences from " << fastaPath << "\n";

    // -------------------------------------------------------------------------
    // Build the DipperSession.
    // -------------------------------------------------------------------------
    DipperSession session;

    // Step 1 (optional) — backbone tree.
    if (withBackbone) {
        std::cerr << "\n=== Mode A: placement onto backbone tree ===\n";
        if (!session.loadBackboneTreeFromFile(backbonePath)) return 1;
    } else {
        std::cerr << "\n=== Mode B: build tree from scratch (no backbone) ===\n";
        std::cerr << "Sequences 0 and 1 will seed the initial tree; "
                  << (seqs.size() - 2) << " queries will be placed.\n";
    }

    // Step 2 — load sequences (backbone + queries, or all sequences).
    session.loadSequences(seqs, names);

    // Step 3 — GPU init.
    session.initializeGPU(params);

    // -------------------------------------------------------------------------
    // Step 4 — place one query sequence at a time using the two-phase API.
    //
    // findNextPlacement()   — DIPPER computes its recommended placement and
    //                          returns a proposal.  The GPU tree is NOT yet
    //                          modified.
    // commitPlacement(a, b) — Insert the query onto the edge between sequential
    //                          labels a and b, updating the GPU tree.
    //                          Pass proposal.splitLabelA / splitLabelB to
    //                          confirm DIPPER's choice, or pass a different pair
    //                          to override and place on a different edge.
    //
    // Sequential labels are pre-order DFS indices starting from 0:
    //   backbone nodes:  0 .. backboneNodeCount()-1
    //   each placement adds:
    //     internalLabel = N   (new internal node)
    //     tipLabel      = N+1 (new tip / query)
    // -------------------------------------------------------------------------
    auto wallStart = std::chrono::high_resolution_clock::now();
    int  placed    = 0;
    int  externalN = (int)session.backboneNodeCount();

    do {
        auto t0 = std::chrono::high_resolution_clock::now();

        // Phase 1: ask DIPPER where it would place the next query.
        PlacementResult proposal = session.findNextPlacement();

        // -----------------------------------------------------------------------
        // *** Integration point ***
        // At this point the GPU tree has NOT been modified yet.
        // You can inspect proposal and update your external node-pointer tree:
        //
        //   proposal.splitLabelA / splitLabelB  — endpoints of the split edge
        //   proposal.internalLabel              — label that the new internal node
        //                                         will receive (regardless of edge choice)
        //   proposal.tipLabel                   — label for the new tip (query)
        //   proposal.tipName                    — sequence name
        //
        // If your external tree topology suggests a different split edge, pass its
        // two endpoint labels to commitPlacement() instead.
        // -----------------------------------------------------------------------

        // Phase 2: commit — here we confirm DIPPER's own recommendation.
        PlacementResult r = session.commitPlacement(
            proposal.splitLabelA,
            proposal.splitLabelB);
        ++placed;

        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        std::cerr << "  query " << placed << " (\"" << r.tipName << "\") placed in "
                  << ms << " ms:\n"
                  << "    split edge  : label " << r.splitLabelA
                  << "  ↔  label " << r.splitLabelB << "\n"
                  << "    new internal: label " << r.internalLabel << "\n"
                  << "    new tip     : label " << r.tipLabel
                  << "  (name: \"" << r.tipName << "\")\n"
                  << "    branch lens : lenA=" << r.lenA
                  << "  lenB=" << r.lenB
                  << "  lenTip=" << r.lenTip << "\n";

        externalN += 2;

        if (!r.more) break;
    } while (true);

    auto wallEnd = std::chrono::high_resolution_clock::now();
    double totalMs =
        std::chrono::duration<double, std::milli>(wallEnd - wallStart).count();
    std::cerr << "Total placement time: " << totalMs << " ms for "
              << placed << " sequence(s).\n";

    // -------------------------------------------------------------------------
    // Output.
    // -------------------------------------------------------------------------
    session.writeTree(outputPath);

    return 0;
}
