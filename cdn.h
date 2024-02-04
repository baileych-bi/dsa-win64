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

#ifndef BIO_CDN_H_
#define BIO_CDN_H_

#include <array>
#include <cassert>
#include <optional>

#include "dna.h"
#include "polymer.h"

namespace bio {

class Cdns;

/** Monomer type where an ASCII char represents a single codon.
  *
  * Unlike DNA and protein, codons have no standard single-letter representation.
  * <pre>To convert an ASCII character codon to nucleotides:
  *   1. subtract 48 from the decimal value (see www.asciitable.com)
  *   2. deconstruct the binary value of resulting byte as follows:
  *      bits 0 and 1 are ignored
  *      bits 2 and 3 encode nucleotide #1
  *      bits 4 and 5 encode nucleotide #2
  *      bits 6 and 7 encode nucleotide #3
  *      According to the following chart:
  *        Binary DNA
  *          00    A 
  *          01    C 
  *          10    T 
  *          11    G
  * Example:
  *   ASCII codon ';' has a decimal value of 59
  *   59 - 48 = 11
  *   11 in binary is 00001011
  *                   ^^        ignored
  *                     ^^      nucleotide #1 is 00 = A
  *                       ^^    nucleotide #2 is 10 = T
  *                         ^^  nucleotide #3 is 11 = G
  *   so ';' is ATG, the start codon.</pre>
  */
class Cdn {
    friend class Cdns;
public:
    static const Cdn AAA;
    static const Cdn AAC;
    static const Cdn AAG;
    static const Cdn AAT;
    static const Cdn ACA;
    static const Cdn ACC;
    static const Cdn ACG;
    static const Cdn ACT;
    static const Cdn AGA;
    static const Cdn AGC;
    static const Cdn AGG;
    static const Cdn AGT;
    static const Cdn ATA;
    static const Cdn ATC;
    static const Cdn ATG;
    static const Cdn ATT;
    static const Cdn CAA;
    static const Cdn CAC;
    static const Cdn CAG;
    static const Cdn CAT;
    static const Cdn CCA;
    static const Cdn CCC;
    static const Cdn CCG;
    static const Cdn CCT;
    static const Cdn CGA;
    static const Cdn CGC;
    static const Cdn CGG;
    static const Cdn CGT;
    static const Cdn CTA;
    static const Cdn CTC;
    static const Cdn CTG;
    static const Cdn CTT;
    static const Cdn GAA;
    static const Cdn GAC;
    static const Cdn GAG;
    static const Cdn GAT;
    static const Cdn GCA;
    static const Cdn GCC;
    static const Cdn GCG;
    static const Cdn GCT;
    static const Cdn GGA;
    static const Cdn GGC;
    static const Cdn GGG;
    static const Cdn GGT;
    static const Cdn GTA;
    static const Cdn GTC;
    static const Cdn GTG;
    static const Cdn GTT;
    static const Cdn TAA;
    static const Cdn TAC;
    static const Cdn TAG;
    static const Cdn TAT;
    static const Cdn TCA;
    static const Cdn TCC;
    static const Cdn TCG;
    static const Cdn TCT;
    static const Cdn TGA;
    static const Cdn TGC;
    static const Cdn TGG;
    static const Cdn TGT;
    static const Cdn TTA;
    static const Cdn TTC;
    static const Cdn TTG;
    static const Cdn TTT;

    Cdn() = default;
    Cdn(const Cdn &) = default;
    Cdn(Nt a, Nt b, Nt c);

    Nt p1() const { return LUT[(v - BIAS) >> 4 & 0b11]; }
    Nt p2() const { return LUT[(v - BIAS) >> 2 & 0b11]; }
    Nt p3() const { return LUT[(v - BIAS) >> 0 & 0b11]; }

    Nt at(size_t i) const { assert(i < 3); return LUT[(v - BIAS) >> (4 - 2*i) & 0b11]; }

    operator std::array<Nt, 3>() const;
    Nts to_nts() const;

    bool operator==(const Cdn &n) const { return v == n.v; }
    bool operator!=(const Cdn &n) const { return v != n.v; }
    bool operator< (const Cdn &n) const { return v <  n.v; }
    bool operator<=(const Cdn &n) const { return v <= n.v; }
    bool operator> (const Cdn &n) const { return v >  n.v; }
    bool operator>=(const Cdn &n) const { return v >= n.v; }

    /** Implicit conversion to char using single ASCII character code explained in class description */
    operator char() const { return v; }

    /** Returns c if c represents a valid codon or 0 otherwise */
    static char normalize_char(char c) { return (AAA.v <= c && c <= GGG.v) ? c : 0; }

    static std::optional<Cdn> from_char(char c) {
        std::optional<Cdn> oc; 
        if (normalize_char(c)) oc = Cdn(c);
        return oc;
    }

    /** Get unique numerical value in range [0..63) */
    size_t index() const { return static_cast<size_t>(v-BIAS); }

    static const std::string valid_chars;

    constexpr static const char BIAS = 0x30;

private:
    static const Nt LUT[4];

    constexpr Cdn(char c) : v(c) {}
    char v = 0x30;
};

class Cdns : public Polymer<Cdn> {
public:
    static const Cdns all;
    static const Cdns all_coding;

    Cdns() = default;

    explicit Cdns(const Cdns &) = default;
    Cdns &operator=(const Cdns &) = default;

    explicit Cdns(Cdns &&cns) { swap_buffers(cns); }
    Cdns &operator=(Cdns &&cns) { swap_buffers(cns); return *this; }

    Cdns(Polymer &&p) : Polymer(std::move(p)) {}

    Cdns(const Nts &dna);
    Cdns &operator=(const Nts &dna);

    Cdns(Nts &&dna);
    Cdns &operator=(Nts &&dna);

    explicit Cdns(const char *s) : Polymer(s) {}
    Cdns &operator=(const char *s);

    explicit Cdns(const std::string &s) : Polymer(s) {}
    Cdns &operator=(const std::string &s);

    Nts to_nts() const;
};

};

template<>
struct std::hash<bio::Cdn>
{
    std::size_t operator()(const bio::Cdn &cdn) const noexcept { 
        return static_cast<size_t>(static_cast<char>(cdn));
    }
};

template<>
struct std::hash<bio::Cdns> {
    std::size_t operator()(const bio::Cdns & cdns) const noexcept {
        return std::hash<std::string_view>{}(cdns.as_string_view());
    }
};


#endif