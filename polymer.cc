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

#include <cassert>
#include <ostream>

#include "polymer.h"

namespace bio {

PolymerBase::PolymerBase(size_t capacity) {
    buf_ = alloc().allocate(capacity, capacity_);
}

PolymerBase::PolymerBase(const PolymerBase &p) {
    assert (&p != this);
    buf_ = alloc().allocate(p.size(), capacity_);
    memcpy(buf_, p.buf_+p.lo_, p.size());
    lo_ = 0;
    hi_ = p.size();
}

PolymerBase::PolymerBase(const char *begin, const char *end) {
    hi_ = end - begin;
    buf_ = alloc().allocate(hi_, capacity_);
    memcpy(buf_, begin, hi_);
}

PolymerBase::~PolymerBase() {
    alloc().deallocate(buf_, capacity_);
}

void
PolymerBase::resize(size_t n, char c) {
    if (n < size()) {
        hi_ -= size() - n;
    } else {
        size_t diff = n - size();
        if (capacity_ - lo_ < n) reserve(n + lo_);
        memset(buf_ + hi_, c, diff);
        hi_ += diff;
    }
}

void
PolymerBase::reserve(size_t n) {
    if (n <= capacity_) return;
    size_t allocd = 0;
    char *tmp = alloc().allocate(n, allocd);
    memcpy(tmp, buf_+lo_, size());
    alloc().deallocate(buf_, 0);
    capacity_ = allocd;
    buf_ = tmp;
    hi_ -= lo_;
    lo_ = 0;
}

void
PolymerBase::shrink_to_fit() {
    if (empty()) {
        alloc().deallocate(buf_, capacity_);
        buf_ = nullptr;
        lo_ = hi_ = capacity_ = 0;
    } else {
        char *tmp = alloc().allocate(size(), capacity_);
        std::memcpy(tmp, buf_, size());
        alloc().deallocate(buf_, capacity_);
        buf_ = tmp;
    }
}

void
PolymerBase::exo(size_t l, size_t r) {
    assert(l+r <= size()); 
    lo_ += l;
    hi_ -= r;
    std::memset(buf_ + hi_, 0, r);
}

void
PolymerBase::pack() {
    PolymerBase tmp(*this);
    swap_buffers(tmp);
}

std::string_view 
PolymerBase::as_string_view() const { 
    return std::string_view(c_data(), size()); 
}

std::ostream &
operator<<(std::ostream &os, const PolymerBase &p) {
    os << p.as_string_view();
    return os;
}

const char *
getline(const char *begin, const char *end, std::string &s) {
    s.clear();
    for (; begin != end; ++begin) {
        const char c = *begin;
        if (c == '\n') return begin + 1;
        s.push_back(c);
    }
    return end;
}

const char *
skipline(const char *begin, const char *end, char delim) {
    for (; begin != end; ++begin) {
        if (*begin == delim) return ++begin;
    }
    return begin;
}

std::istream &
skipline(std::istream &is) {
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
            is.setstate(std::ios::eofbit);
            return is;
        default:
            break;
        }
    }
}

};
