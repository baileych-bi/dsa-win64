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

#ifndef BIO_NTS_H_
#define BIO_NTS_H_

#include <ostream>
#include <string>

#include "polymer.h"

namespace bio {

class Nts;

class Nt {
    friend class Nts;
    friend class Polymer<Nt>;
public:
    static const Nt A;
    static const Nt C;
    static const Nt G;
    static const Nt T;
    static const Nt N;

    Nt() = default;
    Nt(const Nt &) = default;

    bool operator==(const Nt &n) const { return v == n.v; }
    bool operator!=(const Nt &n) const { return v != n.v; }
    bool operator< (const Nt &n) const { return v <  n.v; }
    bool operator<=(const Nt &n) const { return v <= n.v; }
    bool operator> (const Nt &n) const { return v >  n.v; }
    bool operator>=(const Nt &n) const { return v >= n.v; }

    operator char() const { return v; }

    Nt operator ~() const { return Nt(clut[v & 0b1111]); }
    Nt complement() const { return ~*this;               }

    //returns uppercase ATGC for valid or 0 for invalid
    inline static char normalize_char(char c) {
        switch (c) {
            case 'A': case 'a': return 'A';
            case 'T': case 't': return 'T';
            case 'G': case 'g': return 'G';
            case 'C': case 'c': return 'C';
            case 'N': case 'n': return 'N';
            default: return 0;
        };
    };

    static const std::string valid_chars;
    static const char *clut; //complement lookup table: index using (v & 0b1111) to get complement


    static const size_t indexes[8];

    inline size_t index() const { return indexes[static_cast<size_t>((v & 0b1111) >> 1)]; }

private:
    constexpr Nt(char c) : v(c) {}
    char v = 'A';
};

class Cdns;

class Nts : public Polymer<Nt> {
public:
    Nts() = default;
    explicit Nts(size_t capacity) : Polymer(capacity) {}

    Nts(const Nts &) = default;
    Nts &operator=(const Nts &) = default;

    Nts(Nts &&dna) { swap_buffers(dna); }
    Nts &operator=(Nts &&dna) { swap_buffers(dna); return *this; }

    Nts(const Cdns &);
    Nts &operator=(const Cdns &);

    Nts(Polymer &&p) : Polymer(std::move(p)) {}

    Nts(const char *s) : Polymer(s) {}
    Nts(const std::string &s) : Polymer(s) {}

    Nts &complement();
    Nts &reverse_complement();
};

}; //namespace bio

template<>
struct std::hash<bio::Nt> {
    std::size_t operator()(const bio::Nt & nt) const noexcept {
        return static_cast< size_t >( static_cast< char >( nt ) );
    }
};

template<>
struct std::hash<bio::Nts> {
    std::size_t operator()(const bio::Nts & nts) const noexcept {
        return std::hash<std::string_view>{}( nts.as_string_view() );
    }
};

#endif