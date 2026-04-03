# DIPPER Library API — Integration Guide

This guide explains how to embed DIPPER's GPU-accelerated phylogenetic placement
engine as a C++ library inside another tool.

## Table of Contents

1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Building the Examples](#building-the-examples)
4. [Integrating DIPPER Into Your Project](#integrating-dipper-into-your-project)
5. [API Reference](#api-reference)
   - [DipperSession](#dippersession)
   - [PlacementResult](#placementresult)
6. [Workflows](#workflows)
   - [With Backbone Tree](#with-backbone-tree)
   - [Without Backbone Tree](#without-backbone-tree)
7. [Placement APIs](#placement-apis)
   - [One-Phase (Simple)](#one-phase-simple)
   - [Two-Phase (Synchronized)](#two-phase-synchronized)
8. [Node Labeling](#node-labeling)
9. [Complete Code Snippet](#complete-code-snippet)

---

## Overview

`DipperSession` (`src/dipper_api.hpp`) is the single entry point for library
use. It encapsulates all GPU state and exposes a clean step-by-step interface
for:

1. Loading an optional backbone tree (Newick).
2. Loading aligned sequences (FASTA / in-memory strings).
3. Uploading data structures to the GPU.
4. Placing query sequences one at a time, returning structured results after
   each placement.

Two placement modes are available:

| Mode | Backbone? | Behaviour |
|------|-----------|-----------|
| **Backbone placement** | Yes | Query sequences not present in the backbone are placed onto it. |
| **De-novo (from scratch)** | No | First two sequences seed a two-node tree; all remaining sequences are placed in order. |

---

## Requirements

| Dependency | Minimum version |
|------------|----------------|
| CUDA toolkit | 11.0 |
| CMake | 3.10 |
| C++ compiler | GCC 9+ (C++17) |
| Boost | 1.70 (`program_options`) |
| TBB (Intel oneAPI Threading Building Blocks) | 2021 |
| ZLIB | any recent version |

The DIPPER source tree itself has no additional install step — you compile its
source files directly into your executable.

---

## Building the Examples

The examples live in this directory and are built with their own `CMakeLists.txt`
**independently of the main DIPPER project**.

```bash
cd examples
mkdir build && cd build

# CUDA build (required for the GPU-accelerated API)
cmake -DUSE_CUDA=ON ..

# Optional: if your DIPPER source tree is not at ../src, specify it:
# cmake -DUSE_CUDA=ON -DDIPPER_SRC_DIR=/path/to/dipper/src ..

make -j$(nproc)
```

This produces a `dipper_example_msa` binary in the build directory.

### Running the example

```bash
# With a backbone tree:
./dipper_example_msa --backbone backbone.nwk all_sequences.fa output.nwk

# Without a backbone tree (build from scratch):
./dipper_example_msa sequences.fa output.nwk
```

---

## Integrating DIPPER Into Your Project

DIPPER is header-and-source based — there is no pre-built library to link
against. Add its source files to your own build system.

### With CMake

Copy the `add_dipper_example` macro pattern from `examples/CMakeLists.txt` or
adapt it:

```cmake
# In your CMakeLists.txt

set(DIPPER_SRC /path/to/dipper/src)

set(DIPPER_SRCS
    ${DIPPER_SRC}/placement.cu
    ${DIPPER_SRC}/mash.cu
    ${DIPPER_SRC}/twoBitCompressor.cpp
    ${DIPPER_SRC}/fourBitCompressor.cpp
    ${DIPPER_SRC}/matrix_reader.cu
    ${DIPPER_SRC}/MSA.cu
    ${DIPPER_SRC}/tree.cpp
    ${DIPPER_SRC}/neighborJoining.cu
    ${DIPPER_SRC}/placement_close_k.cu
    ${DIPPER_SRC}/divide_and_conquer/mash.cu
    ${DIPPER_SRC}/divide_and_conquer/mash.cpp
    ${DIPPER_SRC}/divide_and_conquer/msa.cu
    ${DIPPER_SRC}/divide_and_conquer/placement_close_k.cpp
    ${DIPPER_SRC}/divide_and_conquer/placement_close_k.cu
)

add_executable(my_tool my_tool.cpp ${DIPPER_SRCS})
set_target_properties(my_tool PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CXX_STANDARD 17)
target_include_directories(my_tool PRIVATE ${DIPPER_SRC})
target_link_libraries(my_tool PRIVATE ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} TBB::tbb)
```

### In your source file

```cpp
#include "dipper_api.hpp"   // only DIPPER header needed
```

---

## API Reference

### DipperSession

```cpp
class DipperSession {
public:
    // ---- Step 1: backbone tree (optional) -----------------------------------
    bool loadBackboneTree(const std::string& newick);
    bool loadBackboneTreeFromFile(const std::string& path);

    // ---- Step 2: sequences --------------------------------------------------
    void loadSequences(const std::vector<std::string>& seqs,
                       const std::vector<std::string>& names);

    // ---- Step 3: GPU initialisation -----------------------------------------
    void initializeGPU(MashPlacement::Param& params);

    // ---- Step 4a: one-phase placement ---------------------------------------
    PlacementResult placeNextSequence();
    void            placeAllSequences();

    // ---- Step 4b: two-phase placement (find, then commit) -------------------
    PlacementResult findNextPlacement();
    PlacementResult commitPlacement(int labelA, int labelB);

    // ---- Output -------------------------------------------------------------
    void writeTree(const std::string& outputPath);
    void writeTree(std::ofstream& ofs);

    // ---- Accessors ----------------------------------------------------------
    size_t backboneSize()      const;   // number of backbone leaf sequences
    size_t backboneNodeCount() const;   // backbone leaves + internal nodes
    size_t numSequences()      const;   // total sequences (backbone + queries)
    bool   hasBackbone()       const;
    int    nextQueryIndex()    const;
    int    nodeLabel(int dipperIdx) const;          // DIPPER idx → sequential label
    const std::vector<std::string>& sortedNames() const;
};
```

### PlacementResult

Returned by `placeNextSequence()`, `findNextPlacement()`, and `commitPlacement()`.

```cpp
struct PlacementResult {
    bool more;          // true if more queries remain to be placed

    // Split edge — the edge that was bisected to insert the new sequence.
    // Before placement:  splitNodeA ──[original len]── splitNodeB
    // After  placement:  splitNodeA ──[lenA]── internal ──[lenB]── splitNodeB
    //                                              └──[lenTip]── tip (query)
    int splitNodeA;     // DIPPER internal index of one endpoint
    int splitNodeB;     // DIPPER internal index of other endpoint
    int splitLabelA;    // sequential label of splitNodeA
    int splitLabelB;    // sequential label of splitNodeB

    // Newly inserted nodes.
    int internalNodeIdx;  // DIPPER internal index of the new internal node
    int tipNodeIdx;       // DIPPER internal index of the new tip (= query seq)
    int internalLabel;    // sequential label of the new internal node
    int tipLabel;         // sequential label of the new tip

    std::string tipName;  // sequence name of the placed query

    // Branch lengths (all ≥ 0).
    double lenA;    // splitNodeA  → internal
    double lenB;    // splitNodeB  → internal
    double lenTip;  // internal    → tip
};
```

---

## Workflows

### With Backbone Tree

```
loadBackboneTree()   ← Newick string or file
loadSequences()      ← FASTA: backbone leaves + query sequences
initializeGPU()
loop:
    findNextPlacement() / placeNextSequence()
    commitPlacement()   (two-phase only)
writeTree()
```

The backbone tree defines the starting topology. Sequences whose names match
backbone leaves are mapped to their existing tree positions; sequences not in
the backbone are treated as queries and placed one at a time.

### Without Backbone Tree

```
loadSequences()    ← FASTA: all sequences (first two seed the tree)
initializeGPU()
loop:
    findNextPlacement() / placeNextSequence()
    commitPlacement()   (two-phase only)
writeTree()
```

The first two sequences in the file seed an initial two-node tree. All
remaining sequences become queries.

---

## Placement APIs

### One-Phase (Simple)

```cpp
do {
    PlacementResult r = session.placeNextSequence();
    // GPU tree is updated immediately.
    // Use r.splitLabelA/B, r.internalLabel, r.tipLabel, r.lenA/B/Tip …
} while (r.more);
```

### Two-Phase (Synchronized)

Use this when your tool maintains its own node-pointer tree alongside DIPPER's
GPU tree and needs to inspect or override the placement before it is applied.

```cpp
do {
    // Phase 1 — DIPPER finds the best edge but does NOT update the GPU tree.
    PlacementResult proposal = session.findNextPlacement();

    // proposal.splitLabelA / splitLabelB  identify the recommended edge.
    // proposal.internalLabel / tipLabel   are the labels the two new nodes
    //                                     will receive regardless of edge choice.
    // proposal.tipName                    is the sequence name.

    // --- Your external tree update goes here --------------------------------
    // e.g. find the edge (proposal.splitLabelA ↔ proposal.splitLabelB) in
    // your tree, split it, and attach internal + tip nodes.
    // ------------------------------------------------------------------------

    // Phase 2 — commit the placement onto the GPU tree.
    // Pass the same labels to confirm DIPPER's choice, or different labels
    // to place on a different edge (override).
    PlacementResult r = session.commitPlacement(
        proposal.splitLabelA,   // ← change these to override the edge
        proposal.splitLabelB);

    // r reflects the edge that was actually used (lengths may differ if overridden).
} while (r.more);
```

**Constraint**: `commitPlacement()` must be called before the next
`findNextPlacement()`. Mixing one-phase and two-phase calls within the same
session is also supported as long as no placement is left uncommitted.

---

## Node Labeling

All backbone nodes are labeled with a sequential pre-order DFS index starting
from 0:

```
backbone root → label 0
  child 1     → label 1
    leaf A    → label 2
    leaf B    → label 3
  leaf C      → label 4
  …
```

`session.backboneNodeCount()` returns the total number of backbone nodes
(leaves + internal), which equals the next available label `N` before any
queries are placed.

Each `placeNextSequence()` / `commitPlacement()` call adds exactly two nodes:

| New node | Label | Meaning |
|----------|-------|---------|
| Internal | `N`   | The bisection point inserted on the split edge |
| Tip      | `N+1` | The newly placed query sequence |

After `k` placements the tree has `backboneNodeCount() + 2k` nodes total.

Use `session.nodeLabel(dipperIdx)` to convert any DIPPER internal index to its
sequential label. Use `session.sortedNames()[dipperIdx]` to get the sequence
name for a given DIPPER index.

---

## Complete Code Snippet

```cpp
#include "dipper_api.hpp"
#include <vector>
#include <string>

void runPlacement(
    const std::string& backboneNewick,          // "" = no backbone
    const std::vector<std::string>& sequences,
    const std::vector<std::string>& names,
    const std::string& outputPath)
{
    MashPlacement::Param params(21, 1000, 0, 2, "m", "t");
    params.range     = {0, -1};
    params.isProtein = false;

    DipperSession session;

    if (!backboneNewick.empty())
        session.loadBackboneTree(backboneNewick);

    session.loadSequences(sequences, names);
    session.initializeGPU(params);

    // backboneNodeCount() is the N to start external node numbering from.
    int N = (int)session.backboneNodeCount();

    do {
        // Phase 1: get DIPPER's recommendation (GPU tree unchanged).
        PlacementResult proposal = session.findNextPlacement();

        // === Update your external tree here using proposal fields ===
        // e.g.:
        //   externalTree.splitEdge(proposal.splitLabelA, proposal.splitLabelB,
        //                          proposal.internalLabel, proposal.tipLabel,
        //                          proposal.lenA, proposal.lenB, proposal.lenTip,
        //                          proposal.tipName);
        N += 2;

        // Phase 2: commit (confirm DIPPER's edge, or pass different labels
        // to override).
        PlacementResult r = session.commitPlacement(
            proposal.splitLabelA,
            proposal.splitLabelB);

        if (!r.more) break;
    } while (true);

    session.writeTree(outputPath);
}
```

See [`msa_placement_example.cpp`](msa_placement_example.cpp) for a complete,
runnable example including FASTA parsing and command-line argument handling.
