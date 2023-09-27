#ifndef BIO_POLYMER_H_
#define BIO_POLYMER_H_

#include "simdalloc.h"

#include <cassert>
#include <iostream>
#include <istream>
#include <iterator>
#include <string>
#include <string_view>

namespace bio {

/**
  * Base class for all polymers (dna, codons, protein)
  *
  * "Polymers" are dynamic arrays of types punnable with char.
  */
class PolymerBase {
    using alloc = ccb::simd_allocator<char, ccb::Register::YMM>;

public:
    PolymerBase() noexcept = default;
    /** Connstruct PolymerBase with 0 size and count capacity.
      *
      * @param count initial capacity
      */
    explicit PolymerBase(size_t count);
    PolymerBase(const PolymerBase &);
    PolymerBase &operator=(const PolymerBase &)=default;

    ~PolymerBase();

    inline size_t size()     const noexcept { return hi_ - lo_; }
    inline size_t capacity() const noexcept { return capacity_; }

    //inline size_t capacity() const { return capacity_ - hi_; }

    void reserve(size_t capacity);

    inline void clear() { hi_ = lo_ = 0; }
    inline bool empty() const { return hi_ == lo_; }

    void shrink_to_fit();

    const char * c_str() const { return buf_ + lo_; } ///< Const access to underlying null-terminated string.
    const char *  data() const { return buf_ + lo_; } ///< Const access to underlying null-terminated string.
    const char *c_data() const { return buf_ + lo_; } ///< Const access to underlying null-terminated string.

          char &front()       { return *(buf_ + lo_); }
    const char &front() const { return *(buf_ + lo_); }

          char &back ()       { return *(buf_ + (hi_ - 1)); }
    const char &back () const { return *(buf_ + (hi_ - 1)); }

    bool operator==(const PolymerBase &rhs) const noexcept { return std::equal(buf_ + lo_, buf_ + hi_, rhs.buf_ + rhs.lo_, rhs.buf_ + rhs.hi_); }

    /** Exonuclease/exoprotease
      *
      * Trim monomers from the left and right. This method never
      * triggers re-allocation of the internal buffer.
      *
      * @param left number of monomers to trim from the left
      * @param right number of monomers to trim from the right
      */
    void exo(size_t left, size_t right);
    void pack();

    std::string_view as_string_view() const;

    friend std::ostream &operator<<(std::ostream &, const PolymerBase &);

protected:
    friend void swap(PolymerBase &a, PolymerBase &b) {
        a.swap_buffers(b);
    }

    explicit PolymerBase(PolymerBase &&p) noexcept { swap_buffers(p); }
    PolymerBase &operator=(PolymerBase &&p) {
        if (this != &p) swap_buffers(p);
        return *this;
    }

    void resize(size_t n, char c);

    void swap_buffers(PolymerBase &p) noexcept {
        std::swap(      lo_, p.lo_);
        std::swap(      hi_, p.hi_);
        std::swap(capacity_, p.capacity_);
        std::swap(     buf_, p.buf_);
    }

    inline void push_back(char c) {
        if (hi_ == capacity_) reserve(capacity_+32);
        buf_[hi_++] = c;
    }

    inline char pop_back() {
        char c = buf_[--hi_];
        buf_[hi_] = 0;
        return c;
    }

    char *data() { return buf_ + lo_; }

    PolymerBase(const char *begin, const char *end);

    size_t lo_ = 0, hi_ = 0, capacity_ = 0;
    char  *buf_ = nullptr;
};

/**
  * Polymer is a dynamic array of type Monomer
  *
  * Many specializations of Polymer support efficient interconversion using
  * move construction and move assignment by re-using the internal buffer in
  * the moved-from object.
  *
  * @tparam Monomer a 1-byte trivially copyable type implicity convertable to a printable ASCII char.
  *         Every class Monomer must have a static method with signature
  *         <pre>static char Monomer::normalize_char(char c);</pre>
  *         that returns a normalized representation of c (e.g. capital ATGC for nucleotides) or 0 for
  *         an invalid/non-convertable character.
  */
template<typename Monomer>
class Polymer : public PolymerBase {
    static_assert (sizeof(Monomer) == sizeof(char), 
        "Polymer: type 'Monomer' must be interconvertible with type 'char' via type punning.");
    static_assert (alignof(Monomer) == alignof(char), 
        "Polymer: type 'Monomer' must be interconvertible with type 'char' via type punning.");
public:
    Polymer() noexcept = default;
    Polymer(const Polymer &) = default;
    Polymer &operator=(const Polymer &) = default;
    Polymer(Polymer &&p) noexcept : PolymerBase(std::move(p)) {}
    Polymer &operator=(Polymer &&) noexcept = default;

    Polymer(const char *);
    Polymer(const std::string &s);

    explicit Polymer(size_t capacity) : PolymerBase(capacity) {}

    void resize(size_t n, const Monomer &m=Monomer());

    /** Append c if c can be converted to a valid Monomer @return true if c was valid, false otherwise */
    inline bool    push_back(char c)    { c=Monomer::normalize_char(c); if (!c) return false; PolymerBase::push_back(c); return true; }
    inline void    push_back(Monomer n) { PolymerBase::push_back(n); }
    inline Monomer pop_back()           { return Monomer(PolymerBase::pop_back()); }

    inline       Monomer &operator[](size_t i);
    inline const Monomer &operator[](size_t i) const;

          Monomer &front()       { return *reinterpret_cast<      Monomer *>(buf_ + lo_); }
    const Monomer &front() const { return *reinterpret_cast<const Monomer *>(buf_ + lo_); }

          Monomer &back ()       { return *reinterpret_cast<      Monomer *>(buf_ + (hi_ - 1)); }
    const Monomer &back () const { return *reinterpret_cast<const Monomer *>(buf_ + (hi_ - 1)); }

    Polymer &operator+=(const Polymer &);

    Polymer subclone(size_t pos, size_t len=std::string::npos);

    template<typename T, typename difference_type>
    struct Forward {
        inline T* inc(T* p, difference_type d)  const { return p+d; }
        inline T* dec(T* p, difference_type d)  const { return p-d; }
        inline T* get(T* p)                     const { return p  ; }
        inline difference_type dist(T *p, T *q) const {return p-q; } 
    };

    template<typename T, typename difference_type>
    struct Reverse {
        inline T* inc(T* p, difference_type d)  const { return p-d; }
        inline T* dec(T* p, difference_type d)  const { return p+d; }
        inline T* get(T* p)                     const { return p-1; }
        inline difference_type dist(T *p, T *q) const {return q-p; } 
    };

    template<bool is_mut, bool forward=true>
    class iter {
        friend class Polymer;
    public:
        using iterator_category = std::contiguous_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value     = typename std::conditional<is_mut, Monomer  , const Monomer  >::type;
        using value_type = value;
        using pointer   = typename std::conditional<is_mut, Monomer *, const Monomer *>::type;
        using reference = typename std::conditional<is_mut, Monomer &, const Monomer &>::type;

        using pun_type  = typename std::conditional<is_mut, char    *, const char    *>::type;

        using advance   = typename std::conditional<forward,
                                                    Forward<value, difference_type>,
                                                    Reverse<value, difference_type>>::type;

        iter() = default;
        iter(const iter &) = default;
        iter(iter<!is_mut,  forward> ii) : p(ii.p) {}
        iter(iter< is_mut, !forward> ii) : p(ii.p) {}

        iter<is_mut, !forward> base() const { return iter<is_mut, !forward>(p); }

        reference operator* () { return *advance().get(p); }
        pointer   operator->() { return  advance().get(p); }
        reference operator[](difference_type n) { return *advance().inc(p, n); }

        iter &operator++(                 ) { p=advance().inc(p, 1); return *this; }
        iter &operator+=(difference_type d) { p=advance().inc(p, d); return *this; }
        iter  operator+ (difference_type d) { return iter(advance().inc(p, d));    }
        iter  operator++(int) { iter tmp = *this; ++(*this); return tmp; }
        friend 
        iter  operator+(const iter &ii, difference_type d) { return ii+d; }

        iter &operator--(                 ) { p=advance().dec(p, 1); return *this; }
        iter &operator-=(difference_type d) { p=advance().dec(p, d); return *this; }
        iter  operator- (difference_type d) { return iter(advance().dec(p, d));    }
        iter  operator--(int) { iter tmp = *this; --(*this); return tmp; }
        friend 
        iter  operator-(const iter &ii, difference_type d) { return ii-d; }

        difference_type operator-(const iter &ii) const { return advance().dist(p, ii.p); }

        friend bool operator==(const iter &a, const iter &b) { return a.p == b.p; };
        friend bool operator!=(const iter &a, const iter &b) { return a.p != b.p; };
        friend bool operator< (const iter &a, const iter &b) { return a.p <  b.p; };
        friend bool operator<=(const iter &a, const iter &b) { return a.p <= b.p; };
        friend bool operator> (const iter &a, const iter &b) { return a.p >  b.p; };
        friend bool operator>=(const iter &a, const iter &b) { return a.p >= b.p; };

    protected:
        explicit iter(pun_type s) : p(reinterpret_cast<pointer>(s)) {}
        explicit iter(pointer  q) : p(q) {}
        pointer p = nullptr;
    };

    using iterator               = iter<true , true >; ///< Random access iterator
    using const_iterator         = iter<false, true >; ///< Random access const iterator
    using reverse_iterator       = iter<true , false>; ///< Random access reverse iterator
    using const_reverse_iterator = iter<false, false>; ///< Random access const reverse iterator

          iterator  begin()       { return       iterator(&buf_[lo_]); }
          iterator  end  ()       { return       iterator(&buf_[hi_]); }
    const_iterator  begin() const { return const_iterator(&buf_[lo_]); }
    const_iterator  end  () const { return const_iterator(&buf_[hi_]); }
    const_iterator cbegin() const { return const_iterator(&buf_[lo_]); }
    const_iterator cend  () const { return const_iterator(&buf_[hi_]); }

    reverse_iterator        rbegin()       { return       reverse_iterator(&buf_[hi_]); }
    reverse_iterator        rend  ()       { return       reverse_iterator(&buf_[lo_]); }
    const_reverse_iterator  rbegin() const { return const_reverse_iterator(&buf_[hi_]); }
    const_reverse_iterator  rend  () const { return const_reverse_iterator(&buf_[lo_]); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(&buf_[hi_]); }
    const_reverse_iterator crend  () const { return const_reverse_iterator(&buf_[lo_]); }

private:
    Polymer(const char *begin, const char *end) : PolymerBase(begin, end) {}
};

template<typename Monomer>
Polymer<Monomer>::Polymer(const char *s) {
    clear();
    for (; *s; ++s) {
        char c = Monomer::normalize_char(*s);
        if (c) PolymerBase::push_back(c);
    }
}

template<typename Monomer>
Polymer<Monomer>::Polymer(const std::string &s) {
    clear();
    reserve(s.size());
    for (char c : s) {
        c = Monomer::normalize_char(c);
        if(c) PolymerBase::push_back(c);
    }
}

template<typename Monomer>
void
Polymer<Monomer>::resize(size_t n, const Monomer &m) {
    PolymerBase::resize(n, static_cast<char>(m));
}

template<typename Monomer>
Monomer &
Polymer<Monomer>::operator[](size_t i) {
    return reinterpret_cast<Monomer &>(buf_[i + lo_]); 
}

template<typename Monomer>
const Monomer &
Polymer<Monomer>::operator[](size_t i) const { 
    return reinterpret_cast<const Monomer &>(buf_[i + lo_]);
}

template<typename Monomer>
Polymer<Monomer> &
Polymer<Monomer>::operator+=(const Polymer<Monomer> &p) {
    if (this != &p) {
        size_t pos = size();
        resize(size() + p.size());
        std::copy(p.begin(), p.end(), begin() + pos);
    } else {
        resize(2*size());
        std::copy(begin(), begin()+size()/2, begin()+size()/2);
    }
    return *this;
}

template<typename Monomer>
Polymer<Monomer>
Polymer<Monomer>::subclone(size_t pos, size_t len) {
    assert (pos <= size());
    Polymer p;
    len = std::min(len, size() - pos);
    p.resize(len);
    std::memcpy(&(p[0]), &((*this)[pos]), len);
    return p;
}

/** Implement getline for polymer types that handles \n, \r, and \r\n newlines.
  * Thanks to Stackoverflow queston #6089231!
  */
template<typename Monomer>
std::istream &
getline(std::istream &is, Polymer<Monomer> &p, size_t &stripped) {
    p.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n') sb->sbumpc();
            return is;
        case std::streambuf::traits_type::eof():
            if(p.empty()) is.setstate(std::ios::eofbit);
            return is;
        default:
            stripped += !p.push_back(static_cast<char>(c));
        }
    }
}

const char *
getline(const char *begin, const char *end, std::string &s);

template<typename Monomer>
const char *
getline(const char *begin, const char *end, Polymer<Monomer> &p, size_t &stripped) {
    p.clear();
    stripped=0;

    for (; begin != end; ++begin) {
        const char c = *begin;
        if (c == '\n') return begin + 1;
        stripped += !p.push_back(static_cast<char>(c));
    }
    return end;
}

const char *
skipline(const char *begin, const char *end, char delim);

/** Ignore a line entirely. */
std::istream &
skipline(std::istream &is);

};
#endif