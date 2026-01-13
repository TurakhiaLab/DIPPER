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

        uint64_t fourBitVal = 0;
        uint64_t shift = 0;
        for (auto j=start; j<end; j++) {
            switch(seq[j+lowerLimit]) {
            case 'A':
                fourBitVal = 0; // DNA/RNA/Protein sequences
                break;
            case 'C':
                fourBitVal = 1; // DNA/RNA/Protein sequences
                break;
            case 'G':
                fourBitVal = 2; // DNA/RNA/Protein sequences
                break;
            case 'T':
                fourBitVal = 3; // DNA/Protein sequences
                break;
            case 'U':
                fourBitVal = 3; // RNA sequences
                break;
            case 'D':
                fourBitVal = 4; // Protein sequences
                break;
            case 'E':
                fourBitVal = 5; // Protein sequences
                break;
            case 'F':
                fourBitVal = 6; // Protein sequences
                break;
            case 'H':
                fourBitVal = 7; // Protein sequences
                break;
            case 'I':
                fourBitVal = 8; // Protein sequences
                break;
            case 'K':
                fourBitVal = 9; // Protein sequences
                break;
            case 'L':
                fourBitVal = 10; // Protein sequences
                break;
            case 'M':
                fourBitVal = 11; // Protein sequences
                break;  
            case 'N':
                fourBitVal = 12; // Protein sequences
                break;
            case 'P':
                fourBitVal = 13; // Protein sequences
                break;
            case 'Q':
                fourBitVal = 14; // Protein sequences
                break;
            case 'R':
                fourBitVal = 15; // Protein sequences
                break;
            case 'S':
                fourBitVal = 16; // Protein sequences
                break;
            case 'V':
                fourBitVal = 17; // Protein sequences
                break;
            case 'W':
                fourBitVal = 18; // Protein sequences
                break;
            case 'Y':
                fourBitVal = 19; // Protein sequences
                break;
            default:
                fourBitVal = isProtein?20:4; // Unknown character
                break;
            }

            compressedSeq[i] |= (fourBitVal << shift);
            if (isProtein) shift += 8;
            else shift += 4;
        }
    }
}