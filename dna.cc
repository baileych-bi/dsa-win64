
#include <immintrin.h>

#include <iostream>

#include "dna.h"
#include "cdn.h"

namespace bio {

const Nt Nt::A = Nt('A');
const Nt Nt::C = Nt('C');
const Nt Nt::G = Nt('G');
const Nt Nt::T = Nt('T');
const Nt Nt::N = Nt('N');

/* ASCII low nibbles    >> 1
 * A: 0b0001 = 1         0    
 * C: 0b0011 = 3         1
 * G: 0b0111 = 7         3
 * T: 0b0100 = 4         2
 * N: 0b1110 = 14        7
 */
const char *Nt::clut = "-T-GA--C------N";
const size_t Nt::indexes[8] = {0, 1, 2, 3, 0, 0, 0, 4};
/* In-place complement the nucleotide sequence in [dna, dna+len) */
void
mm256_complement_dna(char *dna, size_t len) {
    constexpr const size_t chunk = sizeof (__m256i);
    //because Nts contains only A, C, G, T, and N characters
    //we can shuffle based on their low nibbles to get the
    //complement
    __m256i clutv = _mm256_setr_epi8(
        0x80,  'T', 0x80,  'G',  'A', 0x80, 0x80,  'C',
        0x80, 0x80, 0x00, 0x80, 0x80, 0x80,  'N', 0x80,
        0x80,  'T', 0x80,  'G',  'A', 0x80, 0x80,  'C',
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80,  'N', 0x80
    );

    //equivalent lookup table for scalar portion
    size_t i=0;
    for (; i+chunk < len; i += chunk) {
        __m256i seq = _mm256_loadu_si256(reinterpret_cast<__m256i *>(dna+i));
                seq = _mm256_shuffle_epi8(clutv, seq);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dna+i), seq);
    }
    for (; i < len; ++i) *(dna+i) = Nt::clut[static_cast<size_t>(*(dna+i)) & 0x0F];
}

/* In-place reverse-complement the nucleotide sequence in [dna, dna+len) */
void
mm256_reverse_complement_dna(char *dna, size_t len) {
    constexpr const size_t chunk = sizeof(__m256i);
    //because Nts contains only A, C, G, and T characters
    //we can shuffle based on their low nibbles to get the
    //complement
    __m256i clutv = _mm256_setr_epi8(
        0x80,  'T', 0x80,  'G',  'A', 0x80, 0x80,  'C',
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80,  'N', 0x80,
        0x80,  'T', 0x80,  'G',  'A', 0x80, 0x80,  'C',
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80,  'N', 0x80
    );

    __m256i rlutv = _mm256_setr_epi8(
        0x0F, 0x0E, 0x0D, 0x0C, 0x0B, 0x0A, 0x09, 0x08,
        0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00,
        0x0F, 0x0E, 0x0D, 0x0C, 0x0B, 0x0A, 0x09, 0x08,
        0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01, 0x00
    );

    size_t i=0;
    for (; i+2*chunk<=len; i += chunk, len -= chunk) {
        __m256i lseq = _mm256_loadu_si256(reinterpret_cast<__m256i *>(dna+i        ));
                lseq = _mm256_shuffle_epi8(clutv, lseq);
                lseq = _mm256_shuffle_epi8(lseq, rlutv);
                lseq = _mm256_permute4x64_epi64(lseq, 0b01'00'11'10);
        __m256i rseq = _mm256_loadu_si256(reinterpret_cast<__m256i *>(dna+len-chunk));
                rseq = _mm256_shuffle_epi8(clutv, rseq);
                rseq = _mm256_shuffle_epi8(rseq, rlutv);
                rseq = _mm256_permute4x64_epi64(rseq, 0b01'00'11'10);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dna+i          ), rseq);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dna+len-chunk), lseq);
    }
    for (; i+chunk<=len; i+= chunk) {
        __m256i rseq = _mm256_loadu_si256(reinterpret_cast<__m256i *>(dna+len-chunk));
                rseq = _mm256_shuffle_epi8(clutv, rseq);
                rseq = _mm256_shuffle_epi8(rseq, rlutv);
                rseq = _mm256_permute4x64_epi64(rseq, 0b01'00'11'10);
        std::memcpy(dna+i+chunk, dna+i, len-i-chunk);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(dna+i), rseq);
    }
    for (; i+1<=len; ++i, --len) {
        char l = dna[i    ];
        char r = dna[len-1];
        dna[i    ] = Nt::clut[static_cast<size_t>(r) & 0x0F];
        dna[len-1] = Nt::clut[static_cast<size_t>(l) & 0x0F];
    }
}

const std::string
Nt::valid_chars = "ACGTN";

Nts &
Nts::complement() {
    mm256_complement_dna(data(), size());
    return *this;
}

Nts &
Nts::reverse_complement() {
    mm256_reverse_complement_dna(data(), size());
    return *this;
}

Nts::Nts(const Cdns &cdns) : Nts() {
    *this = cdns;
}

Nts &
Nts::operator=(const Cdns &cdns) {
    reserve(3 * cdns.size());
    clear();
    std::array<Nt, 3> nnn;
    for (Cdn c : cdns) {
        nnn = c;
        push_back(nnn[0]);
        push_back(nnn[1]);
        push_back(nnn[2]);
    }

    return *this;
}

}; //namespace bio