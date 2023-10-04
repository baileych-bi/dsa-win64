
#include <immintrin.h>

#include "cdn.h"

namespace bio {

const Cdn Cdn::ACC = Cdn('5');
const Cdn Cdn::ATG = Cdn(';');
const Cdn Cdn::AAG = Cdn('3');
const Cdn Cdn::AAA = Cdn('0');
const Cdn Cdn::ATC = Cdn('9');
const Cdn Cdn::AAC = Cdn('1');
const Cdn Cdn::ATA = Cdn('8');
const Cdn Cdn::AGG = Cdn('?');
const Cdn Cdn::CCT = Cdn('F');
const Cdn Cdn::CTC = Cdn('I');
const Cdn Cdn::AGC = Cdn('=');
const Cdn Cdn::ACA = Cdn('4');
const Cdn Cdn::AGA = Cdn('<');
const Cdn Cdn::CAT = Cdn('B');
const Cdn Cdn::AAT = Cdn('2');
const Cdn Cdn::ATT = Cdn(':');
const Cdn Cdn::CTG = Cdn('K');
const Cdn Cdn::CTA = Cdn('H');
const Cdn Cdn::ACT = Cdn('6');
const Cdn Cdn::CAC = Cdn('A');
const Cdn Cdn::ACG = Cdn('7');
const Cdn Cdn::CAA = Cdn('@');
const Cdn Cdn::AGT = Cdn('>');
const Cdn Cdn::CCA = Cdn('D');
const Cdn Cdn::CCG = Cdn('G');
const Cdn Cdn::CCC = Cdn('E');
const Cdn Cdn::TAT = Cdn('R');
const Cdn Cdn::GGT = Cdn('n');
const Cdn Cdn::TGT = Cdn('^');
const Cdn Cdn::CGA = Cdn('L');
const Cdn Cdn::CAG = Cdn('C');
const Cdn Cdn::CGC = Cdn('M');
const Cdn Cdn::GAT = Cdn('b');
const Cdn Cdn::CGG = Cdn('O');
const Cdn Cdn::CTT = Cdn('J');
const Cdn Cdn::TGC = Cdn(']');
const Cdn Cdn::GGG = Cdn('o');
const Cdn Cdn::TAG = Cdn('S');
const Cdn Cdn::GGA = Cdn('l');
const Cdn Cdn::TAA = Cdn('P');
const Cdn Cdn::GGC = Cdn('m');
const Cdn Cdn::TAC = Cdn('Q');
const Cdn Cdn::TTC = Cdn('Y');
const Cdn Cdn::TCG = Cdn('W');
const Cdn Cdn::TTA = Cdn('X');
const Cdn Cdn::TTG = Cdn('[');
const Cdn Cdn::CGT = Cdn('N');
const Cdn Cdn::TTT = Cdn('Z');
const Cdn Cdn::TCA = Cdn('T');
const Cdn Cdn::GCA = Cdn('d');
const Cdn Cdn::GTA = Cdn('h');
const Cdn Cdn::GCC = Cdn('e');
const Cdn Cdn::GTC = Cdn('i');
const Cdn Cdn::GCG = Cdn('g');
const Cdn Cdn::GTG = Cdn('k');
const Cdn Cdn::GAG = Cdn('c');
const Cdn Cdn::GTT = Cdn('j');
const Cdn Cdn::GCT = Cdn('f');
const Cdn Cdn::TGA = Cdn('\\');
const Cdn Cdn::GAC = Cdn('a');
const Cdn Cdn::TCC = Cdn('U');
const Cdn Cdn::TGG = Cdn('_');
const Cdn Cdn::GAA = Cdn('`');
const Cdn Cdn::TCT = Cdn('V');

const std::string Cdn::valid_chars = "0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmno";
const Cdns Cdns::all(Cdn::valid_chars);
const Cdns Cdns::all_coding("0123456789:;<=>?@ABCDEFGHIJKLMNOQRTUVWXYZ[]^_`abcdefghijklmno");

const Nt Cdn::LUT[4] = {Nt::A, Nt::C, Nt::T, Nt::G};

Cdn::Cdn(Nt a, Nt b, Nt c) {
    v = 0;
    v |= (static_cast<char>(a) & 0b110) << 3;
    v |= (static_cast<char>(b) & 0b110) << 1;
    v |= (static_cast<char>(c) & 0b110) >> 1;
    v += BIAS;
}

Cdn::operator std::array<Nt, 3>() const {
    const char x = v - BIAS;
    std::array<Nt, 3> nnn;
    nnn[0] = LUT[x >> 4 & 0b11];
    nnn[1] = LUT[x >> 2 & 0b11];
    nnn[2] = LUT[x >> 0 & 0b11];
    return nnn;
}

Nts
Cdn::to_nts() const {
    const char x = v - BIAS;
    Nts nnn(3);
    nnn.push_back(LUT[x >> 4 & 0b11]);
    nnn.push_back(LUT[x >> 2 & 0b11]);
    nnn.push_back(LUT[x >> 0 & 0b11]);
    return nnn;
}


/* Codon packing scheme uses bits 5 and 6 from ASCII
A 0100'0001 -> 0
C 0100'0011 -> 1
G 0100'0111 -> 3
T 0101'0100 -> 2
        ^^
*/
void
mm256_pack_cdns(char *dst, const char *src, size_t len) {
    for (size_t i=0, j=0; i<len; i+=30, j+=10) {
        __m256i cdns = _mm256_setzero_si256();
        
        __m256i nts = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src+i));
                nts = _mm256_and_si256(nts, _mm256_set1_epi8(0b0000'0110));

        __m256i lut1 = _mm256_setr_epi8(
            0x00u, 0x03u, 0x06u, 0x09u, 0x0Cu, 0x0Fu, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);

        __m256i lut2 = _mm256_setr_epi8(
            0x01u, 0x04u, 0x07u, 0x0Au, 0x0Du, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);

        __m256i lut3 = _mm256_setr_epi8(
            0x02u, 0x05u, 0x08u, 0x0Bu, 0x0Eu, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);
        
        __m256i p1 = _mm256_shuffle_epi8(nts, lut1);
        __m256i p2 = _mm256_shuffle_epi8(nts, lut2);
        __m256i p3 = _mm256_shuffle_epi8(nts, lut3);

        cdns = _mm256_or_si256(cdns, _mm256_slli_epi16(p1, 3));
        cdns = _mm256_or_si256(cdns, _mm256_slli_epi16(p2, 1));
        cdns = _mm256_or_si256(cdns, _mm256_srli_epi16(p3, 1));

        nts = _mm256_permute4x64_epi64(nts, 0b01'00'11'10);

        lut1 = _mm256_setr_epi8(
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x02u, 0x05u,
            0x08u, 0x0Bu, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);

        lut2 = _mm256_setr_epi8(
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x00u, 0x03u, 0x06u,
            0x09u, 0x0Cu, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);

        lut3 = _mm256_setr_epi8(
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x01u, 0x04u, 0x07u,
            0x0Au, 0x0Du, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u,
            0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u, 0x80u);

        p1 = _mm256_shuffle_epi8(nts, lut1);
        p2 = _mm256_shuffle_epi8(nts, lut2);
        p3 = _mm256_shuffle_epi8(nts, lut3);

        cdns = _mm256_or_si256(cdns, _mm256_slli_epi16(p1, 3));
        cdns = _mm256_or_si256(cdns, _mm256_slli_epi16(p2, 1));
        cdns = _mm256_or_si256(cdns, _mm256_srli_epi16(p3, 1));
 
        __m256i bias = _mm256_set1_epi8(Cdn::BIAS);

        cdns = _mm256_add_epi8(cdns, bias);
        std::memcpy(reinterpret_cast<      void *>(dst + j), 
                    reinterpret_cast<const void *>(&cdns),
                    10);
    }
}

Cdns::Cdns(const Nts &dna) {
    clear();
    resize(dna.size()/3);
    mm256_pack_cdns(data(), dna.c_data(), dna.size());
}

Cdns &
Cdns::operator=(const Nts &dna) {
    clear();
    resize(dna.size()/3);
    mm256_pack_cdns(data(), dna.c_data(), dna.size());
    return *this;
}

Cdns::Cdns(Nts &&dna) {
    swap_buffers(dna);
    mm256_pack_cdns(data(), c_data(), size());
    resize(size()/3);
}

Cdns &
Cdns::operator=(Nts &&dna) {
    swap_buffers(dna);
    mm256_pack_cdns(data(), c_data(), size());
    resize(size()/3);
    return *this;
}

Cdns &
Cdns::operator=(const char *s) {
    clear();
    for ( ; *s; ++s) {
        char c = Cdn::normalize_char(*s);
        if (c) push_back(c);
    }
    return *this;
}

Cdns &
Cdns::operator=(const std::string &s) {
    reserve(s.size());
    return *this = s.c_str();
}

Nts
Cdns::to_nts() const {
    Nts nts(3 * size());
    for (Cdn nnn : *this) {
        nts.push_back(nnn.p1());
        nts.push_back(nnn.p2());
        nts.push_back(nnn.p3());
    }
    return nts;
}

};