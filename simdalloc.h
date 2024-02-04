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

#ifndef CCB_SIMDALLOC_H_
#define CCB_SIMDALLOC_H_

#include <algorithm>
#include <cstring>
#include <string>
#include <cstddef>
#include <iostream>
#include <type_traits>

#ifdef DSA_TARGET_WIN64
#include <malloc.h>
namespace detail {

inline void *
dsa_aligned_alloc(size_t alignment, size_t n) {
    return _aligned_malloc(n, alignment);
}

inline void
dsa_aligned_free(void* p) {
    return _aligned_free(p);
}

}; // namespace detail
#elif defined(DSA_TARGET_LINUX)
#include <cstdlib>
namespace detail {
inline void* dsa_aligned_alloc(size_t alignment, size_t n) {
    return std::aligned_alloc(alignment, n);
}
inline void dsa_aligned_free(void* p) {
    return std::free(p);
}
}; //namespace detail
#endif

namespace ccb {

enum class Register : std::size_t {
    XMM = 16,
    YMM = 32,
    ZMM = 64
};

constexpr const Register REGISTER      = Register::YMM;
constexpr const size_t   REGISTER_SIZE = static_cast<size_t>(REGISTER);

/** 
 * Performs allocation for c++ containers with special guarantees:
 * <ol>
 * <li>memory is aligned to proper boundaries for SIMD instructions (e.g. 32 byte alignment for AVX2/YMM)</li>
 * <li>at least one additional register's worth of bytes is allocated 
 *    beyond what the caller requests (e.g. n+32 bytes for YMM)
 *    unless zero bytes are requested, in which the allocator returns nullptr</li>
 * <li>memory is zero-initialized</li>
 * <li>amount allocated is under-reported by 1 byte, ensuring a null terminator</li>
 * </ol>
 * These guarantees allow many functions to be written purely in terms of SIMD instructions
 * with no need for a scalar loop.
 */
template <typename T, Register Reg> 
class simd_allocator{
public:
    typedef T              value_type;
    typedef std::size_t    size_type;
    typedef std::ptrdiff_t difference_type;
    typedef std::true_type propagate_on_container_move_assignment;

    constexpr  simd_allocator() noexcept = default;
    constexpr  simd_allocator(const simd_allocator &) noexcept = default;
    template<class U> struct rebind { typedef simd_allocator<U, Reg> other; };

    constexpr ~simd_allocator() = default;

    [[nodiscard]] inline constexpr T *allocate(std::size_t n);
     //allocates at least n * sizeof(T) bytes and stores number of bytes allocated in actual
    [[nodiscard]]        constexpr T *allocate(std::size_t n, size_t &actual);

    constexpr void deallocate(T *p, std::size_t n);

};

template <typename T, Register Reg>
constexpr typename simd_allocator<T, Reg>::value_type *
simd_allocator<T, Reg>::allocate(std::size_t n) {
    size_t dummy;
    return allocate(n, dummy);
}

template <typename T, Register Reg>
constexpr typename simd_allocator<T, Reg>::value_type *
simd_allocator<T, Reg>::allocate(std::size_t n, size_t &actual) {
    constexpr const size_t ALN = static_cast<size_t>(Reg);
    actual = 0;
    if (0 == n) return nullptr;
    if ((std::numeric_limits<std::size_t>::max()-ALN) / sizeof(T) < n) throw std::bad_array_new_length();
    actual = (sizeof(T) * n + ALN - 1)/ALN * ALN + ALN;
    void* ptr = detail::dsa_aligned_alloc(ALN, actual);
    if (!ptr) throw std::bad_alloc();
    memset(ptr, 0, actual);
    --actual;
    return reinterpret_cast<value_type *>(ptr);
}

template <typename T, Register Reg>
constexpr void
simd_allocator<T, Reg>::deallocate(value_type *ptr, std::size_t n) {
    (void)n;
    detail::dsa_aligned_free(ptr);
}

}; //namespace ccb

#endif