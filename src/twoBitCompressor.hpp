/**
 * src/twoBitCompressor.hpp
 *
 * 2-bit encoding of DNA/RNA (A=0, C=1, G=2, T/U=3). Used for unaligned sequences
 * prior to MASH sketching. Output: uint64_t array, (seqLen+31)/32 elements.
 * Packs 32 bases per uint64_t; ambiguous bases encoded as 0.
 */

#ifndef TWOBITCOMPRESSOR_HPP
#define TWOBITCOMPRESSOR_HPP

#include <string>
#include <iostream>
#include <tbb/parallel_for.h>
#include "kseq.h"

/** Compress seq[0..seqLen-1] into compressedSeq; caller allocates (seqLen+31)/32 uint64_ts. */
void twoBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq);

#endif