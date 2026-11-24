#ifndef MSA_FILTER_HPP
#define MSA_FILTER_HPP

#include <vector>
#include <string>
#include <cctype>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <tbb/tbb.h>

namespace MsaFilter {

/**
 * Check if a character represents a gap or unknown residue.
 * Consistent with DIPPER's fourBitCompressor encoding.
 */
inline bool isGapChar(char c, bool isProtein) {
    char ch = std::toupper(static_cast<unsigned char>(c));
    if (ch == '-' || ch == '.' || ch == '?') return true;
    if (isProtein) {
        // Standard 20 amino acids: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
        switch (ch) {
            case 'A': case 'C': case 'D': case 'E': case 'F':
            case 'G': case 'H': case 'I': case 'K': case 'L':
            case 'M': case 'N': case 'P': case 'Q': case 'R':
            case 'S': case 'T': case 'V': case 'W': case 'Y':
                return false;
            default:
                return true;
        }
    } else {
        // Standard nucleotides: A, C, G, T, U
        switch (ch) {
            case 'A': case 'C': case 'G': case 'T': case 'U':
                return false;
            default:
                return true;
        }
    }
}

/**
 * Trim aligned MSA columns with gap proportion >= gapThreshold.
 *
 * @param seqs Aligned sequences (modified in-place)
 * @param gapThreshold Gap proportion threshold in [0.0, 1.0] (default: 0.99; >= 1.0 or <= 0.0 disables)
 * @param isProtein True if sequences are amino acids, false for DNA/RNA
 * @return Number of columns trimmed
 */
inline size_t trimGappyColumns(std::vector<std::string>& seqs, double gapThreshold = 0.99, bool isProtein = false) {
    if (seqs.empty() || gapThreshold <= 0.0 || gapThreshold > 1.0) {
        return 0;
    }

    const size_t numSeqs = seqs.size();
    if (numSeqs == 0) return 0;

    const size_t origLen = seqs[0].size();
    if (origLen == 0) return 0;

    // Verify sequences have equal length
    for (size_t i = 1; i < numSeqs; ++i) {
        if (seqs[i].size() != origLen) {
            std::cerr << "Warning: Sequences in MSA have unequal lengths. Skipping gap trimming.\n";
            return 0;
        }
    }

    auto start = std::chrono::high_resolution_clock::now();

    // Parallel count gaps per column using TBB over disjoint column blocks (cache-friendly)
    std::vector<size_t> colGapCounts(origLen, 0);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, origLen), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t row = 0; row < numSeqs; ++row) {
            const std::string& s = seqs[row];
            for (size_t col = range.begin(); col < range.end(); ++col) {
                if (isGapChar(s[col], isProtein)) {
                    colGapCounts[col]++;
                }
            }
        }
    });

    // Select retained columns
    std::vector<size_t> keepCols;
    keepCols.reserve(origLen);
    for (size_t col = 0; col < origLen; ++col) {
        double gapFraction = static_cast<double>(colGapCounts[col]) / static_cast<double>(numSeqs);
        if (gapFraction < gapThreshold) {
            keepCols.push_back(col);
        }
    }

    const size_t trimmedCols = origLen - keepCols.size();
    if (trimmedCols == 0) {
        return 0;
    }

    // Safeguard: do not trim all columns
    if (keepCols.empty()) {
        std::cerr << "Warning: Gap trimming with threshold " << gapThreshold 
                  << " would remove all columns! Trimming skipped.\n";
        return 0;
    }

    // Rebuild sequences with retained columns in parallel
    const size_t newLen = keepCols.size();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, numSeqs), [&](const tbb::blocked_range<size_t>& range) {
        for (size_t row = range.begin(); row < range.end(); ++row) {
            std::string trimmed;
            trimmed.reserve(newLen);
            for (size_t col : keepCols) {
                trimmed.push_back(seqs[row][col]);
            }
            seqs[row] = std::move(trimmed);
        }
    });

    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();

    std::cerr << "Trimmed " << trimmedCols << " / " << origLen 
              << " gappy columns (" << (trimmedCols * 100.0 / origLen) 
              << "%, gap >= " << (gapThreshold * 100.0) << "%). "
              << "Alignment length: " << origLen << " -> " << newLen 
              << " in " << ms << " ms.\n";

    return trimmedCols;
}

} // namespace MsaFilter

#endif // MSA_FILTER_HPP
