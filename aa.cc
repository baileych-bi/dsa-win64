
#include <immintrin.h>
#include "aa.h"

namespace bio {

void
mm256_translate_cdns(char *dst, const char *src, size_t n, const Aa *ttable) {
    constexpr const size_t chunk = sizeof(__m256i);

    __m256i lut0x30 = _mm256_loadu2_m128i(reinterpret_cast<const __m128i *>(ttable+ 0),
                                          reinterpret_cast<const __m128i *>(ttable+ 0));
    __m256i lut0x40 = _mm256_loadu2_m128i(reinterpret_cast<const __m128i *>(ttable+16),
                                          reinterpret_cast<const __m128i *>(ttable+16));
    __m256i lut0x50 = _mm256_loadu2_m128i(reinterpret_cast<const __m128i *>(ttable+32),
                                          reinterpret_cast<const __m128i *>(ttable+32));
    __m256i lut0x60 = _mm256_loadu2_m128i(reinterpret_cast<const __m128i *>(ttable+48),
                                          reinterpret_cast<const __m128i *>(ttable+48));

    __m256i himask = _mm256_set1_epi8(0xF0u);

    for (size_t i=0; i<n; i+= chunk) {
        __m256i cdns   = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src+i));

        __m256i aas    = _mm256_shuffle_epi8(lut0x30, cdns);

        __m256i xlate  = _mm256_shuffle_epi8(lut0x40, cdns);
        __m256i uselut = _mm256_cmpeq_epi8(_mm256_and_si256(cdns, himask), _mm256_set1_epi8(0x040));
                aas    = _mm256_blendv_epi8(aas, xlate, uselut);

                xlate  = _mm256_shuffle_epi8(lut0x50, cdns);
                uselut = _mm256_cmpeq_epi8(_mm256_and_si256(cdns, himask), _mm256_set1_epi8(0x050));
                aas    = _mm256_blendv_epi8(aas, xlate, uselut);

                xlate  = _mm256_shuffle_epi8(lut0x60, cdns);
                uselut = _mm256_cmpeq_epi8(_mm256_and_si256(cdns, himask), _mm256_set1_epi8(0x060));
                aas    = _mm256_blendv_epi8(aas, xlate, uselut);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dst+i), aas);
    }
}

const Aa Aa::A    = Aa('A');
const Aa Aa::C    = Aa('C');
const Aa Aa::D    = Aa('D');
const Aa Aa::E    = Aa('E');
const Aa Aa::F    = Aa('F');
const Aa Aa::G    = Aa('G');
const Aa Aa::H    = Aa('H');
const Aa Aa::I    = Aa('I');
const Aa Aa::K    = Aa('K');
const Aa Aa::L    = Aa('L');
const Aa Aa::M    = Aa('M');
const Aa Aa::N    = Aa('N');
const Aa Aa::P    = Aa('P');
const Aa Aa::Q    = Aa('Q');
const Aa Aa::R    = Aa('R');
const Aa Aa::S    = Aa('S');
const Aa Aa::T    = Aa('T');
const Aa Aa::V    = Aa('V');
const Aa Aa::W    = Aa('W');
const Aa Aa::Y    = Aa('Y');
const Aa Aa::STOP = Aa('*');

const std::string Aa::valid_chars = "*ACDEFGHIKLMNPQRSTVWY";
const Aas Aas::all                = Aa::valid_chars;
const Aas Aas::all_coding         = Aa::valid_chars.substr(1);

const std::array<uint_fast16_t, 48> Aa::indices = {
     0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  1,
     0,  2,  3,  4,  5,  6,  7,  8,
     0,  9, 10, 11, 12,  0, 13, 14,
    15, 16, 17,  0, 18, 19,  0, 20
};

TranslationTable::TranslationTable(const Aas &aas) {
    const size_t n = std::min(table_.size(), aas.size());
    std::copy(aas.begin(), aas.begin()+n, table_.begin());
}

const TranslationTable StandardTranslationTable(Aas("KNNKTTTTIIIMRSSRQHHQPPPPLLLLRRRR*YY*SSSSLFFL*CCWEDDEAAAAVVVVGGGG"));

Aas::Aas(const Cdns &cdns, const TranslationTable &ttable) {
    clear();
    resize(cdns.size());
    mm256_translate_cdns(data(), cdns.c_data(), cdns.size(), ttable.data());
}

void
Aas::set_from_cdns(const Cdns &cdns, const TranslationTable &ttable) {
    clear();
    resize(cdns.size());
    mm256_translate_cdns(data(), cdns.c_data(), cdns.size(), ttable.data());
}

Aas::Aas(Cdns &&cdns,  const TranslationTable &ttable) {
    swap_buffers(cdns);
    mm256_translate_cdns(data(), c_data(), size(), ttable.data());
}

};