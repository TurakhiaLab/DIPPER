#ifndef TWOBITCOMPRESSOR_HPP
#include "twoBitCompressor.hpp"
#endif

void twoBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq) {
    size_t compressedSeqLen = (seqLen+31)/32;

    for (size_t i=0; i < compressedSeqLen; i++) {
        compressedSeq[i] = 0;

        size_t start = 32*i;
        size_t end = std::min(seqLen, start+32);

        uint64_t twoBitVal = 0;
        uint64_t shift = 0;
        for (auto j=start; j<end; j++) {
            switch(seq[j]) {
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
                twoBitVal = 3;
                break;
            default:
                twoBitVal = 0;
                break;
            }

            compressedSeq[i] |= (twoBitVal << shift);
            shift += 2;
        }
    }
}