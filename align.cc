/*
Copyright 2024, The Broad Institute of MIT and Harvard

Original Author: Charles C Bailey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <immintrin.h>

#include <iomanip>
#include <unordered_map>

#include "align.h"
#include "aa.h"
#include "cdn.h"

namespace bio {

Overlap
find_overlapv_256(const char *a, const size_t a_size, const char *b, const size_t b_size, size_t max_mismatches) {
    constexpr const unsigned REGW = sizeof(__m256i);
    std::vector<uint16_t> upper((a_size+1+REGW-1)/REGW*REGW+REGW, 0);
    std::vector<uint16_t> lower((a_size+1+REGW-1)/REGW*REGW+REGW, 0);

    bool in_order = true;
    size_t max_overlap = 0, max_row = 0;
    for (size_t r=0; r<b_size; ++r) {
        __m256i row_char = _mm256_set1_epi8(b[r]);
        std::swap(upper, lower);
        for (size_t c=0; c<a_size; c+=REGW) {
            __m256i col_chars   = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&a[c]));
            __m256i match       = _mm256_cmpeq_epi8(col_chars, row_char);       //0xFF where match
                    match       = _mm256_and_si256(match, _mm256_set1_epi8(1)); //0x01 where match
            __m256i lmatch      = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(match, 0));
            __m256i rmatch      = _mm256_cvtepi8_epi16(_mm256_extracti128_si256(match, 1));
            __m256i scores      = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&lower[c   ]));
                    scores      = _mm256_adds_epu16(scores, lmatch);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(&upper[c+1   ]), scores);
                    scores      = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(&lower[c+16]));
                    scores      = _mm256_adds_epu16(scores, rmatch);
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(&upper[c+1+16]), scores);
        }
        bool new_max = (max_overlap < upper[a_size]) && (r + 1 <= upper[a_size] + max_mismatches);
        max_overlap = max_overlap * (!new_max) + new_max * upper[a_size];
        max_row = max_row * (!new_max) + new_max * r;
    }
    
    for (size_t c=0; c <a_size; ++c) {
        bool new_max = (max_overlap < upper[c+1]) && (c+1 <= upper[c+1] + max_mismatches);
        max_overlap  = max_overlap * (!new_max) + new_max * upper[c+1];
        max_row      = max_row     * (!new_max) + new_max * c;
        in_order     = in_order * !new_max;
    }
    return {max_row+1, max_row+1-max_overlap, in_order};
}

Read
Read::assemble(Read &&fw, Read &&rv, size_t min_overlap_size, size_t max_mismatches) {
    Read rd;

    rv.dna.reverse_complement();

    Overlap ol = find_overlapv_256(fw.dna.c_data(), fw.dna.size(), 
                                   rv.dna.c_data(), rv.dna.size());
    if (ol.overlap < min_overlap_size || ol.mismatches > max_mismatches) return rd;

    std::reverse(rv.qual.begin(), rv.qual.end());

    if (!ol.in_order) {
        std::swap(fw.dna,  rv.dna);
        std::swap(fw.qual, rv.qual);
    }

    for (size_t i=fw.dna.size()-ol.overlap, j=0; j<ol.overlap; ++i, ++j) {
        if (fw.qual[i] < rv.qual[j]) {
            fw.qual[i] = rv.qual[j];
            fw.dna[i]  = rv.dna[j] ;
        }
    }

    rv.dna.exo(ol.overlap, 0);
    fw.dna += rv.dna;
    fw.qual.append(rv.qual.begin()+ol.overlap, rv.qual.end());

    rd.barcode = std::move(fw.barcode);
    rd.barcode += rv.barcode;
    rd.dna  = std::move(fw.dna );
    rd.qual = std::move(fw.qual);

    return rd;
}

static const std::vector<int32_t> blosum62_data = {
 0,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
-4,4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2,
-4,0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,
-4,-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3,
-4,-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2,
-4,-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3,
-4,0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3,
-4,-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2,
-4,-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1,
-4,-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2,
-4,-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1,
-4,-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1,
-4,-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2,
-4,-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3,
-4,-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1,
-4,-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2,
-4,1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2,
-4,0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2,
-4,0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1,
-4,-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2,
-4,-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7
};

static const std::vector<int32_t> cdnsubs_data = {
6,0,0,5,-1,-1,-1,-1,-3,-3,-3,-1,2,0,0,2,1,-1,-1,1,-1,-1,-1,-1,-2,-2,-2,-2,2,2,2,2,-4,-2,-2,-4,0,0,0,0,-2,-3,-3,-2,-4,-3,-3,-3,1,-1,-1,1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
0,7,6,0,0,0,0,0,-3,-3,-3,-2,0,1,1,0,0,1,1,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0,-4,-2,-2,-4,1,1,1,1,-3,-3,-3,-3,-4,-3,-3,-4,0,1,1,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0, 
0,6,7,0,0,0,0,0,-3,-3,-3,-2,0,1,1,0,0,1,1,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0,-4,-2,-2,-4,1,1,1,1,-3,-3,-3,-3,-4,-3,-3,-4,0,1,1,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0, 
5,0,0,6,-1,-1,-1,-1,-3,-3,-3,-1,2,0,0,2,1,-1,-1,1,-1,-1,-1,-1,-2,-2,-2,-2,2,2,2,2,-4,-2,-2,-4,0,0,0,0,-2,-3,-3,-2,-4,-3,-3,-3,1,-1,-1,1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,0,0,-1,6,5,5,5,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,-1,-1,-2,-1,-1,-1,-1,0,0,0,0,0,0,0,0,-2,-2,-2,-2, 
-1,0,0,-1,5,6,5,5,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,-1,-1,-2,-1,-1,-1,-1,0,0,0,0,0,0,0,0,-2,-2,-2,-2, 
-1,0,0,-1,5,5,6,5,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,-1,-1,-2,-1,-1,-1,-1,0,0,0,0,0,0,0,0,-2,-2,-2,-2, 
-1,0,0,-1,5,5,5,6,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,-1,-1,-2,-1,-1,-1,-1,0,0,0,0,0,0,0,0,-2,-2,-2,-2, 
-3,-3,-3,-3,-1,-1,-1,-1,5,4,4,1,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,2,2,2,2,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,2,0,0,2,-4,-1,-1,-3,-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,-4,-4,-4,-4, 
-3,-3,-3,-3,-1,-1,-1,-1,4,5,4,1,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,2,2,2,2,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,2,0,0,2,-4,-1,-1,-3,-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,-4,-4,-4,-4, 
-3,-3,-3,-3,-1,-1,-1,-1,4,4,5,1,-3,-2,-2,-3,-3,-3,-3,-3,-3,-3,-3,-3,2,2,2,2,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,2,0,0,2,-4,-1,-1,-3,-3,-3,-3,-3,-1,-1,-1,-1,3,3,3,3,-4,-4,-4,-4, 
-1,-2,-2,-1,-1,-1,-1,-1,1,1,1,6,-1,-1,-1,-1,0,-2,-2,0,-2,-2,-2,-2,2,2,2,2,-1,-1,-1,-1,-4,-1,-1,-4,-1,-1,-1,-1,2,0,0,2,-4,-1,-1,-1,-2,-3,-3,-2,-1,-1,-1,-1,1,1,1,1,-3,-3,-3,-3, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,6,-1,-1,5,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,5,5,5,5,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,5,4,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,4,4,4,4,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,4,5,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,4,4,4,4,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,5,-1,-1,6,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,5,5,5,5,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
1,0,0,1,-1,-1,-1,-1,-3,-3,-3,0,1,0,0,1,6,0,0,5,-1,-1,-1,-1,-2,-2,-2,-2,1,1,1,1,-4,-1,-1,-4,0,0,0,0,-2,-3,-3,-2,-4,-3,-3,-2,2,0,0,2,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,1,1,-1,-2,-2,-2,-2,-3,-3,-3,-2,0,-1,-1,0,0,9,8,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0,-4,2,2,-4,-1,-1,-1,-1,-3,-1,-1,-3,-4,-3,-3,-2,0,-1,-1,0,-2,-2,-2,-2,-3,-3,-3,-3,-2,-2,-2,-2, 
-1,1,1,-1,-2,-2,-2,-2,-3,-3,-3,-2,0,-1,-1,0,0,8,9,0,-2,-2,-2,-2,-3,-3,-3,-3,0,0,0,0,-4,2,2,-4,-1,-1,-1,-1,-3,-1,-1,-3,-4,-3,-3,-2,0,-1,-1,0,-2,-2,-2,-2,-3,-3,-3,-3,-2,-2,-2,-2, 
1,0,0,1,-1,-1,-1,-1,-3,-3,-3,0,1,0,0,1,5,0,0,6,-1,-1,-1,-1,-2,-2,-2,-2,1,1,1,1,-4,-1,-1,-4,0,0,0,0,-2,-3,-3,-2,-4,-3,-3,-2,2,0,0,2,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,-2,-2,-1,-1,-1,-1,-1,-3,-3,-3,-2,-2,-1,-1,-2,-1,-2,-2,-1,8,7,7,7,-3,-3,-3,-3,-2,-2,-2,-2,-4,-3,-3,-4,-1,-1,-1,-1,-3,-4,-4,-3,-4,-3,-3,-4,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,-2,-2,-1,-1,-1,-1,-1,-3,-3,-3,-2,-2,-1,-1,-2,-1,-2,-2,-1,7,8,7,7,-3,-3,-3,-3,-2,-2,-2,-2,-4,-3,-3,-4,-1,-1,-1,-1,-3,-4,-4,-3,-4,-3,-3,-4,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,-2,-2,-1,-1,-1,-1,-1,-3,-3,-3,-2,-2,-1,-1,-2,-1,-2,-2,-1,7,7,8,7,-3,-3,-3,-3,-2,-2,-2,-2,-4,-3,-3,-4,-1,-1,-1,-1,-3,-4,-4,-3,-4,-3,-3,-4,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,-2,-2,-1,-1,-1,-1,-1,-3,-3,-3,-2,-2,-1,-1,-2,-1,-2,-2,-1,7,7,7,8,-3,-3,-3,-3,-2,-2,-2,-2,-4,-3,-3,-4,-1,-1,-1,-1,-3,-4,-4,-3,-4,-3,-3,-4,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,5,4,4,4,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,4,0,0,4,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,4,5,4,4,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,4,0,0,4,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,4,4,5,4,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,4,0,0,4,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,4,4,4,5,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,4,0,0,4,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,5,-1,-1,5,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,6,5,5,5,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,5,-1,-1,5,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,5,6,5,5,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,5,-1,-1,5,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,5,5,6,5,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-1,5,-1,-1,5,1,0,0,1,-2,-2,-2,-2,-2,-2,-2,-2,5,5,5,6,-4,-2,-2,-4,-1,-1,-1,-1,-2,-3,-3,-2,-4,-3,-3,-3,0,-2,-2,0,-1,-1,-1,-1,-3,-3,-3,-3,-2,-2,-2,-2, 
-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1,-4,-4,0,-4,-4,-4,-4,-4,-4,-4,-4,0,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 
-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-2,-2,-2,-2,-1,2,2,-1,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2,-2,-4,8,7,-4,-2,-2,-2,-2,-1,3,3,-1,-4,-2,-2,2,-2,-3,-3,-2,-2,-2,-2,-2,-1,-1,-1,-1,-3,-3,-3,-3, 
-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-2,-2,-2,-2,-1,2,2,-1,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2,-2,-4,7,8,-4,-2,-2,-2,-2,-1,3,3,-1,-4,-2,-2,2,-2,-3,-3,-2,-2,-2,-2,-2,-1,-1,-1,-1,-3,-3,-3,-3, 
-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,0,-4,-4,1,-4,-4,-4,-4,-4,-4,-4,-4,0,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,4,4,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,5,4,4,4,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,4,4,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,4,5,4,4,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,4,4,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,4,4,5,4,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
0,1,1,0,1,1,1,1,-2,-2,-2,-1,-1,4,4,-1,0,-1,-1,0,-1,-1,-1,-1,-2,-2,-2,-2,-1,-1,-1,-1,-4,-2,-2,-4,4,4,4,5,-2,-2,-2,-2,-4,-1,-1,-3,0,0,0,0,1,1,1,1,-2,-2,-2,-2,0,0,0,0, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,4,4,4,4,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,5,0,0,4,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
-3,-3,-3,-3,-2,-2,-2,-2,0,0,0,0,-3,-2,-2,-3,-3,-1,-1,-3,-4,-4,-4,-4,0,0,0,0,-3,-3,-3,-3,-4,3,3,-4,-2,-2,-2,-2,0,7,6,0,-4,-2,-2,1,-3,-3,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-3,-3,-3,-3, 
-3,-3,-3,-3,-2,-2,-2,-2,0,0,0,0,-3,-2,-2,-3,-3,-1,-1,-3,-4,-4,-4,-4,0,0,0,0,-3,-3,-3,-3,-4,3,3,-4,-2,-2,-2,-2,0,6,7,0,-4,-2,-2,1,-3,-3,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-3,-3,-3,-3, 
-2,-3,-3,-2,-1,-1,-1,-1,2,2,2,2,-2,-2,-2,-2,-2,-3,-3,-2,-3,-3,-3,-3,4,4,4,4,-2,-2,-2,-2,-4,-1,-1,-4,-2,-2,-2,-2,4,0,0,5,-4,-1,-1,-2,-3,-4,-4,-3,-1,-1,-1,-1,1,1,1,1,-4,-4,-4,-4, 
-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,0,-4,-4,0,-4,-4,-4,-4,-4,-4,-4,-4,1,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 
-3,-3,-3,-3,-1,-1,-1,-1,-1,-1,-1,-1,-3,-1,-1,-3,-3,-3,-3,-3,-3,-3,-3,-3,-1,-1,-1,-1,-3,-3,-3,-3,-4,-2,-2,-4,-1,-1,-1,-1,-1,-2,-2,-1,-4,10,9,-2,-4,-3,-3,-4,0,0,0,0,-1,-1,-1,-1,-3,-3,-3,-3, 
-3,-3,-3,-3,-1,-1,-1,-1,-1,-1,-1,-1,-3,-1,-1,-3,-3,-3,-3,-3,-3,-3,-3,-3,-1,-1,-1,-1,-3,-3,-3,-3,-4,-2,-2,-4,-1,-1,-1,-1,-1,-2,-2,-1,-4,9,10,-2,-4,-3,-3,-4,0,0,0,0,-1,-1,-1,-1,-3,-3,-3,-3, 
-3,-4,-4,-3,-2,-2,-2,-2,-3,-3,-3,-1,-3,-3,-3,-3,-2,-2,-2,-2,-4,-4,-4,-4,-2,-2,-2,-2,-3,-3,-3,-3,-4,2,2,-4,-3,-3,-3,-3,-2,1,1,-2,-4,-2,-2,12,-3,-4,-4,-3,-3,-3,-3,-3,-3,-3,-3,-3,-2,-2,-2,-2, 
1,0,0,1,-1,-1,-1,-1,-3,-3,-3,-2,0,0,0,0,2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-3,0,0,0,0,-4,-2,-2,-4,0,0,0,0,-3,-3,-3,-3,-4,-4,-4,-3,6,2,2,5,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,1,1,-1,-1,-1,-1,-1,-3,-3,-3,-3,-2,0,0,-2,0,-1,-1,0,-1,-1,-1,-1,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-4,2,7,6,2,-2,-2,-2,-2,-3,-3,-3,-3,-1,-1,-1,-1, 
-1,1,1,-1,-1,-1,-1,-1,-3,-3,-3,-3,-2,0,0,-2,0,-1,-1,0,-1,-1,-1,-1,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-4,2,6,7,2,-2,-2,-2,-2,-3,-3,-3,-3,-1,-1,-1,-1, 
1,0,0,1,-1,-1,-1,-1,-3,-3,-3,-2,0,0,0,0,2,0,0,2,-1,-1,-1,-1,-3,-3,-3,-3,0,0,0,0,-4,-2,-2,-4,0,0,0,0,-3,-3,-3,-3,-4,-4,-4,-3,5,2,2,6,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,-2,-2, 
-1,-2,-2,-1,0,0,0,0,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,0,0,-3,-1,-2,-2,-1,5,4,4,4,0,0,0,0,0,0,0,0, 
-1,-2,-2,-1,0,0,0,0,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,0,0,-3,-1,-2,-2,-1,4,5,4,4,0,0,0,0,0,0,0,0, 
-1,-2,-2,-1,0,0,0,0,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,0,0,-3,-1,-2,-2,-1,4,4,5,4,0,0,0,0,0,0,0,0, 
-1,-2,-2,-1,0,0,0,0,-1,-1,-1,-1,-1,1,1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,-2,-2,-4,1,1,1,1,-1,-2,-2,-1,-4,0,0,-3,-1,-2,-2,-1,4,4,4,5,0,0,0,0,0,0,0,0, 
-2,-3,-3,-2,0,0,0,0,3,3,3,1,-3,-2,-2,-3,-2,-3,-3,-2,-2,-2,-2,-2,1,1,1,1,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,1,-1,-1,1,-4,-1,-1,-3,-2,-3,-3,-2,0,0,0,0,5,4,4,4,-3,-3,-3,-3, 
-2,-3,-3,-2,0,0,0,0,3,3,3,1,-3,-2,-2,-3,-2,-3,-3,-2,-2,-2,-2,-2,1,1,1,1,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,1,-1,-1,1,-4,-1,-1,-3,-2,-3,-3,-2,0,0,0,0,4,5,4,4,-3,-3,-3,-3, 
-2,-3,-3,-2,0,0,0,0,3,3,3,1,-3,-2,-2,-3,-2,-3,-3,-2,-2,-2,-2,-2,1,1,1,1,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,1,-1,-1,1,-4,-1,-1,-3,-2,-3,-3,-2,0,0,0,0,4,4,5,4,-3,-3,-3,-3, 
-2,-3,-3,-2,0,0,0,0,3,3,3,1,-3,-2,-2,-3,-2,-3,-3,-2,-2,-2,-2,-2,1,1,1,1,-3,-3,-3,-3,-4,-1,-1,-4,-2,-2,-2,-2,1,-1,-1,1,-4,-1,-1,-3,-2,-3,-3,-2,0,0,0,0,4,4,4,5,-3,-3,-3,-3, 
-2,0,0,-2,-2,-2,-2,-2,-4,-4,-4,-3,-2,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-2,-2,-1,-1,-2,0,0,0,0,-3,-3,-3,-3,7,6,6,6, 
-2,0,0,-2,-2,-2,-2,-2,-4,-4,-4,-3,-2,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-2,-2,-1,-1,-2,0,0,0,0,-3,-3,-3,-3,6,7,6,6, 
-2,0,0,-2,-2,-2,-2,-2,-4,-4,-4,-3,-2,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-2,-2,-1,-1,-2,0,0,0,0,-3,-3,-3,-3,6,6,7,6, 
-2,0,0,-2,-2,-2,-2,-2,-4,-4,-4,-3,-2,0,0,-2,-2,-2,-2,-2,-2,-2,-2,-2,-4,-4,-4,-4,-2,-2,-2,-2,-4,-3,-3,-4,0,0,0,0,-4,-3,-3,-4,-4,-3,-3,-2,-2,-1,-1,-2,0,0,0,0,-3,-3,-3,-3,6,6,6,7
};

Matrix<int32_t> print_cdnsubs() {
    Matrix<int32_t> cdnsubs(64, 64);
    for (size_t i=0; i<64; ++i) {
        Cdn cdn1 = *Cdn::from_char(Cdn::valid_chars[i]);
        Aa   aa1 = StandardTranslationTable.translate(cdn1);
        for (size_t j=0; j<64; ++j) {
            Cdn cdn2 = *Cdn::from_char(Cdn::valid_chars[j]);
            Aa aa2 = StandardTranslationTable.translate(cdn2);
            int32_t score = BLOSUM62.elem(aa1.index(), aa2.index());
            score += cdn1 == cdn2;
            cdnsubs.elem(cdn1.index(), cdn2.index()) = score;
            std::cout << score << ", ";
        }
        std::cout << std::endl;
    }
    return cdnsubs;
}

const Matrix<int32_t> BLOSUM62(21, 21, blosum62_data);
const Matrix<int32_t> NTSUBS(4, 4, { 1, -1, -1, -1,
                                    -1,  1, -1, -1,
                                    -1, -1,  1, -1,
                                    -1, -1, -1,  1});
const Matrix<int32_t> CDNSUBS(64, 64, cdnsubs_data);

void
Alignment::clear() {
    score = 0;
    traceback.resize(0, 0);
    aligned_query.clear();
}

void
Read::reverse_complement() {
    resize(size() / 3 * 3);
    dna.reverse_complement();
    std::reverse(qual.begin(), qual.end());
}

Orf::Orf(Read &&rd) {
    umi_group_size = rd.umi_group_size;
    barcode = std::move(rd.barcode);
    cdns = std::move(rd.dna);
    aas.set_from_cdns(cdns);
}

bool
Orf::contains_ptc() const {
    return std::find(aas.begin(), aas.end(), Aa::STOP) != aas.end();
}

/*
void
nw_align_aas(const Aas &q, const Aas &t, const Matrix<int32_t> &match, int32_t gapp, Alignment &result) {
    result.clear();

    thread_local Matrix<Cell> & trace = result.traceback;
    trace.resize(q.size() + 1, t.size() + 1);

    for (int i = 1; i < trace.rows(); ++i ) {
        trace.elem(i, 0).score = 0; //-gapp * i;
        trace.elem(i, 0).move = Cell::Move::GAP_T;
    }
    for (int j = 1; j < trace.cols(); ++j ) {
        trace.elem(0, j).score = 0; //-gapp * j;
        trace.elem(0, j).move = Cell::Move::GAP_Q;
    }

    for ( size_t i = 0; i < q.size(); ++i ) {
        size_t n = q [ i ].index();
        for ( size_t j = 0; j < t.size(); ++j ) {
            size_t m = t [ j ].index();

            Cell cell;
            cell.move = Cell::Move::MATCH;
            cell.score = trace.elem(i, j).score + match.elem(m, n);

            int32_t gappa_score = trace.elem(i + 1, j).score - gapp * (1 + j != t.size());
            if ( gappa_score > cell.score ) {
                cell.score = gappa_score;
                cell.move = Cell::Move::GAP_Q;
            }

            int32_t gappb_score = trace.elem(i, j + 1).score - gapp * (1 + i != q.size());
            if ( gappb_score > cell.score ) {
                cell.score = gappb_score;
                cell.move = Cell::Move::GAP_T;
            }

            trace.elem(i + 1, j + 1) = cell;
        }
    }
    result.score = trace.elem(q.size(), t.size()).score;

    size_t i = q.size(), j = t.size();
    while ( i + j != 0 ) {
        switch ( trace.elem(i, j).move ) {
        case Cell::Move::GAP_Q:
            result.aligned_query.push_back('-');
            --j;
            break;
        case Cell::Move::GAP_T:
            result.aligned_query.push_back(std::tolower(static_cast< char >( q [ i - 1 ] )));
            --i;
            break;
        default:
            result.aligned_query.push_back(std::toupper(static_cast< char >( q [ i - 1 ] )));
            --i;
            --j;
        }
    }
    std::reverse(result.aligned_query.begin(), result.aligned_query.end());
}
*/

std::ostream &
operator<<(std::ostream &os, const Read &rd) {
    os << std::setw(9) << std::left << "barcode"    << " \t" << rd.barcode << std::endl;
    os << std::setw(9) << std::right << "umi count" << " \t" << rd.umi_group_size << std::endl;
    os << std::setw(9) << std::right << "fw dna"    << " \t" << rd.dna << std::endl;
    os << std::setw(9) << std::right << "fw qual"   << " \t" << rd.qual << std::endl; 

    return os;
}

}; //namespace bio

