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

    Nts(Nts &&dna) = default;
    Nts &operator=(Nts &&) = default;

    Nts(const Cdns &);
    Nts &operator=(const Cdns &);

    Nts(Polymer &&p) : Polymer(std::move(p)) {}

    Nts(const char *s) : Polymer(s) {}
    Nts(const std::string &s) : Polymer(s) {}

    Nts &complement();
    Nts &reverse_complement();
};

}; //namespace bio

#endif