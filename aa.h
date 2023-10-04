#ifndef BIO_AA_H_
#define BIO_AA_H_

#include "cdn.h"
#include "polymer.h"

#include <array>
#include <optional>

namespace bio {

class Aa;
class Aas;
class TranslationTable;

/** Monomer type representing amino acids. */
class Aa {
    friend class Aas;
    friend class TranslationTable;
public:
    static const Aa A;
    static const Aa C;
    static const Aa D;
    static const Aa E;
    static const Aa F;
    static const Aa G;
    static const Aa H;
    static const Aa I;
    static const Aa K;
    static const Aa L;
    static const Aa M;
    static const Aa N;
    static const Aa P;
    static const Aa Q;
    static const Aa R;
    static const Aa S;
    static const Aa T;
    static const Aa V;
    static const Aa W;
    static const Aa Y;
    static const Aa STOP;

    Aa() = default; ///< Initialized to stop codon by default.

    /** Implicit conversion to char, capital one-letter amino acid code */
    operator char() const { return v; }

    /** Packed, numeric value of the residue, in the range [0 .. 21) */
    size_t index() const { return indices[static_cast<size_t>(v-'*')]; }; 

    /**
      * Get normalized representation of an amino acid character.
      *
      * For amino acids, this is the capital single letter IUPAC code, or '*' for a stop.
      * Uppercase amino acid codes and '*' will be returned as-is, lower case amino acid codes will 
      * be capitalized, All other characters (invalid amino acids) will be mapped to 0.
      * 
      * @param c an ASCII character
      * @return the normalized character or 0 for an invalid amino acid
      */
    static char 
    normalize_char(char c) { c = std::toupper(c); return (valid_chars.find(c) == std::string::npos) ? 0 : c; }

    /**
      * Attempt char to Aa conversion.
      *
      * @param c an ASCII character
      * @return optional<Aa> containing nullopt if c is not a valid amino acid code.
      */
    static std::optional<Aa> 
    from_char(char c) { std::optional<Aa> aa; if ((c = normalize_char(c))) aa = Aa(c); return aa; }

    static const std::string valid_chars; ///< The list of valid amino acid codes.

private:
    static const std::array<uint_fast16_t, 48> indices;
    constexpr Aa(char c) : v(c) {}
    char v = '*';
};

class TranslationTable {
public:
    TranslationTable() = default;
    TranslationTable(const Aas &);

    Aa translate(Cdn c) const { return table_[c.index()]; }
    const Aa *data() const { return &table_[0]; }

    void map(Cdn c, Aa a) { table_[static_cast<char>(c)] = a; }

private:
    alignas(ccb::REGISTER_SIZE) std::array<Aa, 64> table_;
};

extern const TranslationTable StandardTranslationTable;

class Aas : public Polymer<Aa> {
public:
    static const Aas all;                 ///< A list of all the aas, including *
    static const Aas all_coding;          ///< A list of all the aas, excluding *


    Aas() = default;

    Aas(const Aas &) = default;
    Aas &operator=(const Aas &) = default;

    Aas(Polymer<Aa> &&p) { swap_buffers(p); }

    Aas(Aas &&aas) { swap_buffers(aas); }
    Aas &operator=(Aas &&aas) noexcept { swap_buffers(aas); return *this; }

    Aas(const Cdns &, const TranslationTable &ttable=StandardTranslationTable);
    Aas(Cdns &&, const TranslationTable &ttable=StandardTranslationTable);

    void set_from_cdns(const Cdns &, const TranslationTable &ttable=StandardTranslationTable);

    Aas(Nts &&dna, const TranslationTable &ttable=StandardTranslationTable) : Aas(Cdns(std::move(dna)), ttable) {}

    Aas(const char        *s) : Polymer(s) {}
    Aas(const std::string &s) : Polymer(s) {}
};

};

template<>
struct std::hash<bio::Aa> {
    std::size_t operator()(const bio::Aa &aa) const noexcept {
        return static_cast<size_t>(static_cast<char>(aa));
    }
};

template<>
struct std::hash<bio::Aas> {
    std::size_t operator()(const bio::Aas & aas) const noexcept {
        return std::hash<std::string_view>{}(aas.as_string_view());
    }
};

#endif