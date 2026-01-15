#ifndef FOURBITCOMPRESSOR_HPP
#include "fourBitCompressor.hpp"
#endif

void fourBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq, int lowerLimit, int upperLimit, bool isProtein) {
    
    if (upperLimit > -1) seqLen=upperLimit+1;
    else upperLimit=-1;
    if (lowerLimit > 0) seqLen-=lowerLimit;
    else lowerLimit=0;
    
    size_t compressedSeqLen;
    if (isProtein) compressedSeqLen = (seqLen+7)/8;
    else compressedSeqLen = (seqLen+15)/16;
    
    for (size_t i=0; i < compressedSeqLen; i++) {
        compressedSeq[i] = 0;

        size_t start = isProtein ? 8*i : 16*i;
        size_t end = std::min(seqLen, start + (isProtein ? 8 : 16));

        uint64_t twoBitVal = 0;
        uint64_t shift = 0;
        for (auto j=start; j<end; j++) {
            switch(seq[j+lowerLimit]) {
            case 'A':
                twoBitVal = 0; // DNA/RNA/Protein sequences
                break;
            case 'C':
                twoBitVal = 1; // DNA/RNA/Protein sequences
                break;
            case 'G':
                twoBitVal = 2; // DNA/RNA/Protein sequences
                break;
            case 'T':
                twoBitVal = 3; // DNA/Protein sequences
                break;
            case 'U':
                twoBitVal = 3; // RNA sequences
                break;
            case 'D':
                twoBitVal = 4; // Protein sequences
                break;
            case 'E':
                twoBitVal = 5; // Protein sequences
                break;
            case 'F':
                twoBitVal = 6; // Protein sequences
                break;
            case 'H':
                twoBitVal = 7; // Protein sequences
                break;
            case 'I':
                twoBitVal = 8; // Protein sequences
                break;
            case 'K':
                twoBitVal = 9; // Protein sequences
                break;
            case 'L':
                twoBitVal = 10; // Protein sequences
                break;
            case 'M':
                twoBitVal = 11; // Protein sequences
                break;  
            case 'N':
                twoBitVal = 12; // Protein sequences
                break;
            case 'P':
                twoBitVal = 13; // Protein sequences
                break;
            case 'Q':
                twoBitVal = 14; // Protein sequences
                break;
            case 'R':
                twoBitVal = 15; // Protein sequences
                break;
            case 'S':
                twoBitVal = 16; // Protein sequences
                break;
            case 'V':
                twoBitVal = 17; // Protein sequences
                break;
            case 'W':
                twoBitVal = 18; // Protein sequences
                break;
            case 'Y':
                twoBitVal = 19; // Protein sequences
                break;
            default:
                twoBitVal = 20; // Unknown character
                break;
            }

            compressedSeq[i] |= (twoBitVal << shift);
            if (isProtein) shift += 8;
            else shift += 4;
        }
    }
}