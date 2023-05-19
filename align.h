#ifndef BIO_ALIGN_H_
#define BIO_ALIGN_H_

#include "aa.h"
#include "dna.h"

#include <algorithm>
#include <cassert>
#include <ostream>
#include <optional>
#include <variant>
#include <vector>

namespace bio {

/** Simple 2D array class with 1D backing buffer */
template<typename T>
class Matrix {
public:
    Matrix() = default;
    Matrix(const Matrix &m) = default;
    Matrix(size_t rows, size_t cols, const T &t=T{});

    /** Init from 1D vector.
      * rows * cols must equal init.size()
      */
    Matrix(size_t rows, size_t cols, const std::vector<T> &init);

    Matrix &operator=(const Matrix &m) = default;

    /** Element-wise sum */
    Matrix operator+(const Matrix &m) const;

    /** Reize to fit rows, cols and <em>clear</em> all contents!
      * All elements will be set to t
      */
    void resize(size_t rows, size_t cols, const T &t=T{});

    inline       T &elem(size_t row, size_t col);
    inline const T &elem(size_t row, size_t col) const;
    void fill(const T &);

    size_t rows() const { return rows_; }
    size_t cols() const { return cols_; }

private:
    size_t rows_ = 0, cols_ = 0;
    std::vector<T> buf_;
};

extern const Matrix<int32_t> BLOSUM62; //< The BLOSUM62 matrix for Aa alignments
extern const Matrix<int32_t> NTSUBS;   //< A matrix for Nt alignments, 1 on diagonals and -1 on off-diagonals
extern const Matrix<int32_t> CDNSUBS;  //< A matrix for Cdn alignments, equal to BLOSUM matrix scores +1 for identical codons

Matrix<int32_t> print_cdnsubs();

typedef std::string Qual;

/** Find the length of the longest suffix of a that is also a prefix of b,
  * permitting up to max_mismatches mismatches.
  * Returns length of the longest overlapping region and number of mismatches in that region.
*/
struct Overlap {
    size_t overlap    = 0;     ///< length of overlapping region
    size_t mismatches = 0;     ///< number of mismatches in overlapping region
    /**
      * This is true if sequences overlap 3' to 3' or false if the overlap is 5' to 5'
      */
    bool   in_order   = true;
};

/** Find the length of the longest suffix of a that is also a prefix of b
  * 
  * Vectorized implementation using AVX instructions. Tolerates some number of
  * mismatches in the alignment although for most uses max_mismatches should
  * be left at 0 (the default).
  *
  * @param a sequence a
  * @param a_size number of symbols in a
  * @param b sequence b
  * @param b_size number of symbols in b
  * @param max_mismatches maximum permissible mismatces in overlap region
  * @return Overlap object containing length of the longest overlapping region and number of mismatches in that region.
*/
Overlap
find_overlapv_256(const char *a, 
                  const size_t a_size,
                  const char *b,
                  const size_t b_size,
                  size_t max_mismatches=0);

/** Element of the Needleman-Wunsch traceback matrix. */
struct Cell {
    /** Shows path taken to reach a current cell. */
    enum class Move : int8_t { 
        MATCH = 0, ///< Sequences were matched.
        GAP_Q = 1, ///< A gap was introduced in the query.
        GAP_T = 2  ///< A gap was introduced in the template.
    };
    int32_t score = 0;       ///< The maximum possible alignment score at this location.
    Move move = Move::MATCH; ///< The path taken from the previous Cell.
};

/** The results of a Needleman-Wunsch alignment. */
struct Alignment {
    int32_t score = 0;         ///< The global alignment score.
    Matrix<Cell> traceback;    ///< The traceback matrix.
    std::string aligned_query; ///< A gapped string showing elements in the query that align with those in the template.

    void clear();              ///< Reset for use in another call to nw_align_aas or other function.

    template<typename M>
    std::string build_string(const Polymer<M> &q) const {
        return build_string<M>(q.cbegin(), q.cend());
    }

    template<typename M>
    std::string build_string(typename Polymer<M>::const_iterator q_lo,
                             typename Polymer<M>::const_iterator q_hi) const;
};

template<typename M>
inline char gap_char() { return '-'; }

template<>
inline char gap_char<Cdn>() { return ' '; }

template<typename M>
inline char ins_char(const M &m) { return std::tolower(static_cast<char>(m)); }

template<>
inline char ins_char<Cdn>(const Cdn &c) { return static_cast<char>(c); }

template<typename M>
inline char reg_char(const M &m) { return std::toupper(static_cast<char>(m)); }

template<>
inline char reg_char<Cdn>(const Cdn &c) { return static_cast<char>(c); }

/*
template<typename M>
std::string
Alignment::build_string(typename Polymer<M>::const_iterator q_lo,
                        typename Polymer<M>::const_iterator q_hi,
                        typename Polymer<M>::const_iterator t_lo,
                        typename Polymer<M>::const_iterator t_hi) const {
    const size_t q_size = q_hi - q_lo;
    const size_t t_size = t_hi - t_lo;

    assert (q_size+1 == traceback.rows() && t_size+1 == traceback.cols());

    std::string aligned_query;

    size_t i = q_size, j = t_size;
    while (i + j != 0) {
        switch (traceback.elem(i, j).move) {
            case Cell::Move::GAP_Q:
                aligned_query.push_back(gap_char<M>());
                --j;
            break;
            case Cell::Move::GAP_T:
                aligned_query.push_back(ins_char<M>(q_lo[i-1]));
                --i;
            break;
            default:
                aligned_query.push_back(reg_char<M>(q_lo[i-1]));
                --i;
                --j;
        }
    }
    std::reverse(aligned_query.begin(), aligned_query.end());

    return aligned_query;
}
*/

template<typename M>
std::string
Alignment::build_string(typename Polymer<M>::const_iterator q_lo,
                        typename Polymer<M>::const_iterator q_hi) const {
    const size_t q_size = q_hi - q_lo;
    const size_t t_size = traceback.cols()-1;

    assert (q_size+1 == traceback.rows());

    std::string aligned_query;

    size_t i = q_size, j = t_size;
    while (i + j != 0) {
        switch (traceback.elem(i, j).move) {
            case Cell::Move::GAP_Q:
                aligned_query.push_back(gap_char<M>());
                --j;
            break;
            case Cell::Move::GAP_T:
                aligned_query.push_back(ins_char<M>(q_lo[i-1]));
                --i;
            break;
            default:
                aligned_query.push_back(reg_char<M>(q_lo[i-1]));
                --i;
                --j;
        }
    }
    std::reverse(aligned_query.begin(), aligned_query.end());

    return aligned_query;
}

/** Implementation of Needleman-Wunsch for amino acid sequences.
  *
  * Takes a user-supplied match score matrix and gap penalty. Results are stored
  * in the pass-by-reference Alignment structure. The substitution matrix should 
  * be 21x21 with rows and columns corresponding to the values of Aa::index().
  *
  * @param query the amino acid "query" sequence
  * @param templ the amino acid "template" sequence
  * @param substitution_matrix the substituion scoring matrix (e.g. BLOSUM62)
  * @param gapp the gap extension penalty
  * @param result stores the results of the alignment
  */
void
nw_align_aas(const Aas &query, 
             const Aas &templ, 
             const Matrix<int32_t> &substitution_matrix, 
             int32_t gapp,
             Alignment &result);

/** Generic implementation of Needleman-Wunsch for Polymers of monomer type M.
  *
  * Takes a user-supplied match score matrix and gap penalty. Results are stored
  * in the pass-by-reference Alignment structure. The substitution matrix should 
  * be 4x4 with rows and columns corresponding to the values of Nt::index().
  *
  * @tparam M the type of monomer (Aa, Nt, etc.)
  * @param q_lo the first monomer in the "query" sequence
  * @param q_hi the last + 1 monomer in the "query" sequence
  * @param t_lo the first monomer in the "template" sequence
  * @param t_hi the last + 1 monomer in the "template" sequence
  * @param templ the amino acid "template" sequence
  * @param substitution_matrix the substituion scoring matrix (e.g. BLOSUM62 for M=Aa)
  * @param gapp the gap penalty
  * @param result stores the results of the alignment as a string
  * @param score_only if set to true, only the alignment score will be calculated
  */
template<typename M>
void
nw_align(typename Polymer<M>::const_iterator q_lo,
         typename Polymer<M>::const_iterator q_hi,
         typename Polymer<M>::const_iterator t_lo,
         typename Polymer<M>::const_iterator t_hi,
         const Matrix<int32_t> &match,
         int32_t gapp,
         Alignment &result,
         bool score_only=false) {
    result.clear();

    const size_t q_size = q_hi - q_lo;
    const size_t t_size = t_hi - t_lo;
    thread_local Matrix<Cell> &trace = result.traceback;
    trace.resize(q_size+1, t_size+1);

    for (size_t i = 1; i < trace.rows(); ++i) trace.elem(i, 0).move = Cell::Move::GAP_T; //Cell::Move::GAP_B;
    for (size_t j = 1; j < trace.cols(); ++j) trace.elem(0, j).move = Cell::Move::GAP_Q; //Cell::Move::GAP_A;

    for (size_t i = 0; i < q_size; ++i) {
        size_t n = q_lo[i].index();
        int32_t gapp_a = static_cast<int32_t>(i != q_size-1) * gapp;
        for (size_t j = 0; j < t_size; ++j) {
            size_t m = t_lo[j].index();
            int32_t gapp_b = static_cast<int32_t>(j != t_size-1) * gapp;

            Cell cell;
            cell.move  = Cell::Move::MATCH;
            cell.score = trace.elem(i, j).score + match.elem(m, n);

            int32_t gappa_score = trace.elem(i+1, j).score - gapp_a;
            if (gappa_score > cell.score) {
                cell.score = gappa_score;
                cell.move = Cell::Move::GAP_Q;
            }

            int32_t gappb_score = trace.elem(i, j+1).score - gapp_b;
            if (gappb_score > cell.score) {
                cell.score = gappb_score;
                cell.move  = Cell::Move::GAP_T;
            }

            trace.elem(i+1, j+1) = cell;
        }
    }
    result.score = trace.elem(q_size, t_size).score;

    if (!score_only) result.aligned_query = result.build_string<M>(q_lo, q_hi);
}

/** Generic implementation of Needleman-Wunsch for Polymers of monomer type M.
  *
  * Takes a user-supplied match score matrix and gap penalty. Results are stored
  * in the pass-by-reference Alignment structure. The substitution matrix should 
  * be 4x4 with rows and columns corresponding to the values of Nt::index().
  *
  * @tparam M the type of monomer (Aa, Nt, etc.)
  * @param query the amino acid "query" sequence
  * @param templ the amino acid "template" sequence
  * @param substitution_matrix the substituion scoring matrix (e.g. BLOSUM62 for M=Aa)
  * @param gapp the gap penalty
  * @param result stores the results of the alignment as a string
  * @param score_only if set to true, only the alignment score will be calculated
  */
template<typename M>
void
nw_align(const Polymer<M> &q,
         const Polymer<M> &t,
         const Matrix<int32_t> &match,
         int32_t gapp,
         Alignment &result,
         bool score_only=false) {
    return nw_align<M>(q.cbegin(), q.cend(), t.cbegin(), t.cend(), match, gapp, result, score_only);
}

/** Compute the Needleman-Wunsch score for aligning a sequence to itself. */
template<typename M>
int32_t
nw_self_align_score(const Polymer<M> &query, const Matrix<int32_t> &matrix) {
    int32_t score = 0;
    for (const M &m : query) {
        int32_t index = static_cast<int32_t>(m.index());
        score += matrix.elem(index, index); 
    }
    return score;
}


/** Represents a deep sequencing read data.
  *
  * Could be forward read, reverse read, an assembled read pair, or a consensus
  * sequence of multiple Reads following UMI collapse.
  */
struct Read {
    std::string barcode;       ///< The extracted UMI barcode.
    size_t umi_group_size = 1; ///< The number of reads used to form the consensus.
    Nts  dna;                  ///< The nucleotide sequence.
    Qual qual;                 ///< The fastq quality scores of the nucleotides.

    /** Assemble paired end reads.
      *
      * If assemly fails, the empty() will be true for the returned Read.
      *
      * @param fw the forward read
      * @param rv the reverse read
      * @param min_overlap the minimum acceptible 3' overlap for assembly
      * @param max_mismatches the max allowed mismatches in the 3' overlap region
      * @return a Read obect containing the assembled sequences; empty() will be true if assembly fails
      */
    static Read assemble(
        Read &&fw,
        Read &&rv,
        size_t min_overlap,
        size_t max_mismatches=0);

    bool empty() const { return dna.empty(); }
    size_t size() const { return dna.size(); }
    void resize(size_t n) { dna.resize(n); qual.resize(n); }
    void pop_back() { dna.pop_back(); qual.pop_back(); }
    void reverse_complement();
};

/** A pair of unassembled forward and reverse reads */
struct ReadPair {
    Read fw;
    Read rv;
};

/** Print a read for debugging purposes */
std::ostream &
operator<<(std::ostream &os, const Read &rd);

/** The translation of a Read */
struct Orf {
    Orf() = default;
    Orf(Read &&);

    size_t umi_group_size = 1;
    size_t template_id = 0;
    std::string barcode;
    Cdns cdns;
    Aas  aas;

    bool contains_ptc() const;
};

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T& t)
    : rows_(rows)
    , cols_(cols)
    , buf_(rows*cols, t) {}

template<typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, const std::vector<T> &init)
    : rows_(rows)
    , cols_(cols)
    , buf_(rows*cols) {
        assert(buf_.size() == init.size());
        std::copy(init.begin(), init.end(), buf_.begin());
    }

template<typename T>
void
Matrix<T>::resize(size_t rows, size_t cols, const T& t) {
    rows_ = rows;
    cols_ = cols;
    buf_.clear();
    buf_.resize(rows*cols, t);
}

template<typename T>
T &
Matrix<T>::elem(size_t row, size_t col) {
    return buf_[row*cols_+col];
}

template<typename T>
const T &
Matrix<T>::elem(size_t row, size_t col) const {
    return buf_[row*cols_+col];
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(const Matrix<T> &m) const {
    assert(rows() == m.rows() && cols() == m.cols());
    Matrix<T> sum(*this);
    for (size_t i=0; i<buf_.size(); ++i) sum.buf_[i] += m.buf_[i];
    return sum;
}

template<typename T>
void
Matrix<T>::fill(const T &t) {
    std::fill(buf_.begin(), buf_.end(), t);
}

};

#endif
