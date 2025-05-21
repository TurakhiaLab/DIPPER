#ifndef FOURBITCOMPRESSOR_HPP
#define FOURBITCOMPRESSOR_HPP

#include <string>
#include <iostream>
#include <tbb/parallel_for.h>
#include "kseq.h"

void fourBitCompressor(std::string seq, size_t seqLen, uint64_t* compressedSeq);

#endif