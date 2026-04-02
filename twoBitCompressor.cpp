/**
 * src/twoBitCompressor.cpp
 *
 * 2-bit encoding of DNA/RNA (A=0, C=1, G=2, T/U=3). Packs 32 bases per uint64_t;
 * used for unaligned sequences prior to MASH sketching. Output length
 * (seqLen+31)/32. Ambiguous bases (default) encoded as 0.
 */

#ifndef TWOBITCOMPRESSOR_HPP
#include "twoBitCompressor.hpp"
#endif

/** Compress seq[0..seqLen-1] into compressedSeq; caller allocates (seqLen+31)/32 uint64_ts. */
void twoBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq) {
    size_t compressedSeqLen = (seqLen + 31) / 32;

    /* Process each 64-bit word (32 bases). */
    for (size_t i = 0; i < compressedSeqLen; i++) {
        compressedSeq[i] = 0;
        size_t start = 32 * i;
        size_t end = std::min(seqLen, start + 32);
        uint64_t shift = 0;
        for (auto j = start; j < end; j++) {
            uint64_t twoBitVal;
            switch (seq[j]) {
            case 'A':
                twoBitVal = 0;
                break;
            case 'C':
                twoBitVal = 1;
                break;
            case 'G':
                twoBitVal = 2;
                break;
            case 'T':
            case 'U':
                twoBitVal = 3;
                break;
            default:
                twoBitVal = 0;  /* Ambiguous / unknown: treat as A. */
                break;
            }
            compressedSeq[i] |= (twoBitVal << shift);
            shift += 2;
        }
    }
}