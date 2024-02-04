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

#ifndef CCB_PARALLELISM_H_
#define CCB_PARALLELISM_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <unordered_map>
#include <thread>
#include <vector>

/** 
  * Multithreaded implementations of some algorithms.
  */
namespace bio {

namespace impl {

template<typename InputIt, typename UnaryFunction>
void
for_each(InputIt first, InputIt last, UnaryFunction f) {
    std::for_each(first, last, f);
}

template<typename InputIt, typename OutputT, typename UnaryOperation>
void
transform(InputIt first, InputIt last, std::vector<OutputT> &result, UnaryOperation unary_op) {
    result.clear();
    result.reserve(std::distance(first, last));
    std::transform(first, last, std::back_inserter(result), unary_op);
}

template<typename InputIt, typename OutputT, typename TransformFilter, typename Log>
void
transform_filter(InputIt first, InputIt last, std::vector<OutputT> &result, TransformFilter tf, Log &log) {
    result.clear();
    result.reserve(std::distance(first, last));
    auto bi = std::back_inserter(result);
    for (; first != last; ++first) {
        auto optional_result = tf(*first, log);
        if (optional_result.has_value()) *bi++ = std::move(*optional_result);
    }
}

template<typename InputIt, typename OutputT, typename TransformFilterLog>
void
transform_filter(InputIt first, InputIt last, std::vector<OutputT> &result, TransformFilterLog &tfl) {
    result.clear();
    result.reserve(std::distance(first, last));
    auto bi = std::back_inserter(result);
    for (; first != last; ++first) {
        auto optional_result = tfl(*first);
        if (optional_result.has_value()) *bi++ = std::move(*optional_result);
    }
}

template<typename InputIt, typename Reduce, typename OutputT>
void
reduce(InputIt first, InputIt last, Reduce f, OutputT &out) {
    out = f(first, last);
}

}; //namespace impl

template<typename InputIt>
auto
modal_element(InputIt first, InputIt last)->InputIt {
    if (first == last) return last;

    using key_type = decltype(*first);
    struct counter {
        InputIt iter;
        size_t  count = 0;
    };

    std::unordered_map<key_type, counter> counts;

    using value_type = decltype(*counts.first());

    for (; first != last; ++first) {
        value_type &val = counts[*first];
        val.iter = first;
        val.count += 1;
    }

    struct cmp {
        inline bool operator()(const value_type &a, const value_type &b) const {
            return a.second.count < b.second.count;
        }
    };

    return std::max_element(counts.begin(), counts.end(), cmp())->second.iter;
}


template<typename InputIt, typename UnaryFunction>
void
parallel_for_each(InputIt first, InputIt last, UnaryFunction f) {
    const typename std::iterator_traits<InputIt>::difference_type n = std::distance(first, last);
    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t batch = n / thread_count;

    if (batch == 0) return impl::for_each(first, last, f);

    std::vector<std::thread> threads(thread_count-1);

    InputIt lo = first, hi = lo;
    size_t i = 0;
    for (; i != thread_count-1; ++i) {
        std::advance(hi, batch);
        threads[i] = std::thread(impl::for_each<InputIt, UnaryFunction>, lo, hi, f);
        std::advance(lo, batch);
        hi = lo;
    }
    impl::for_each(lo, last, f);

    for (auto &th : threads) th.join();
}

template<typename InputIt, typename OutputIt, typename UnaryOperation>
void
parallel_transform(InputIt first, InputIt last, OutputIt out, UnaryOperation unary_op) {
    typedef typename std::iterator_traits<InputIt>::value_type InputT;
    typedef decltype(unary_op(*first)) OutputT;

    const size_t n = std::distance(first, last);
    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t batch = n / thread_count;

    if (batch == 0) {
        std::transform(first, last, out, unary_op);
        return;
    }

    std::vector<std::thread> threads(thread_count-1);
    std::vector<std::vector<OutputT>> fragments(thread_count);

    InputIt lo = first, hi = lo;

    size_t i = 0;
    for (; i != thread_count-1; ++i) {
        std::advance(hi, batch);
        threads[i] = std::thread(impl::transform<InputIt, OutputT, UnaryOperation>, 
            lo, hi, std::ref(fragments[i]), unary_op);
        std::advance(lo, batch);
        hi = lo;
    }
    impl::transform(lo, last, fragments[i], unary_op);

    for (auto &th : threads) th.join();
    for (auto &frag : fragments) {
        out = std::copy(
            std::make_move_iterator(frag.begin()),
            std::make_move_iterator(frag.end()),
            out);
    }
}

/**
  * Multithreaded transform over [first, last) using binary functor TransformFilter
  * TransformFilter should return std::optional<T> where T is the desired output type
  * std::nullopt outputs from TransformFilter are discarded
  * 
  * @tparam InputIt the input iterator type
  * @tparam OutputIt the output iterator type (threads std::copy their results here)
  * @tparam TransformFilter a functor that takes InputIt::value_type and Log and returns std::optional<OutputIt::value_type>
  * @tparam Log a structure to capture error information, must define operator + for use in std::accumulate
  * @param first the first iterator in the range
  * @param last the one-past-the end iterator in the range
  * @param out the destination for results
  * @param tf callable that takes InputIt::value_type and a Log instance, returns std::optional<OutputIt::value_type>
  * @return OutputIt out after recieving transformed values
  */
template<typename InputIt, typename OutputIt, typename TransformFilter, typename Log>
OutputIt
parallel_transform_filter(InputIt first, InputIt last, OutputIt out, TransformFilter tf, Log &log) {
    typedef typename std::iterator_traits<InputIt>::value_type InputT;
    typedef typename decltype(tf(*first, log))::value_type OutputT;

    const size_t n = std::distance(first, last);
    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t batch = n / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    std::vector<std::vector<OutputT>> fragments(thread_count);
    std::vector<Log> logs(thread_count);

    InputIt batch_last = first;

    size_t i = 0;
    for (; i != thread_count-1; ++i) {
        std::advance(batch_last, batch);
        threads[i] = std::thread(impl::transform_filter<InputIt, OutputT, TransformFilter, Log>, 
            first, batch_last, std::ref(fragments[i]), tf, std::ref(logs[i]));
        std::advance(first, batch);
        batch_last = first;
    }
    impl::transform_filter(first, last, fragments[i], tf, logs[i]);

    for (auto &th : threads) th.join();
    for (auto &frag : fragments) {
        out = std::copy(
            std::make_move_iterator(frag.begin()),
            std::make_move_iterator(frag.end()),
            out);
    }
    log = std::accumulate(logs.begin(), logs.end(), log);

    return out;
}

template<typename InputIt, typename OutputIt, typename TransformFilterLog>
OutputIt
parallel_transform_filter(InputIt first, InputIt last, OutputIt out, TransformFilterLog &tfl) {
    typedef typename std::iterator_traits<InputIt>::value_type InputT;
    typedef typename decltype(tfl(*first))::value_type OutputT;

    const size_t n = std::distance(first, last);
    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t batch = n / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    std::vector<std::vector<OutputT>> fragments(thread_count);
    std::vector<TransformFilterLog> logs(thread_count);

    InputIt batch_last = first;

    size_t i = 0;
    for (; i != thread_count-1; ++i) {
        std::advance(batch_last, batch);
        threads[i] = std::thread(impl::transform_filter<InputIt, OutputT, TransformFilterLog>, 
            first, batch_last, std::ref(fragments[i]), std::ref(logs[i]));
        std::advance(first, batch);
        batch_last = first;
    }
    impl::transform_filter(first, last, fragments[i], logs[i]);

    for (auto &th : threads) th.join();
    for (auto &frag : fragments) {
        out = std::copy(
            std::make_move_iterator(frag.begin()),
            std::make_move_iterator(frag.end()),
            out);
    }
    tfl = std::accumulate(logs.begin(), logs.end(), tfl);

    return out;
}

/**
  * Multithreaded reduction over [first, last) using binary function Reduces
  * @tparam InputIt the input iterator type
  * @tparam Reduce a binary function that operates on two <em>iterators</em> representing a range of values
  * @param first the first iterator in the range
  * @param last the one-past-the end iterator in the range
  * @param f a callable (funnctor, lambda, etc.) that takes two InputIt
  * @return the return type of f
  */
template<typename InputIt, typename Reduce>
auto
parallel_reduce(InputIt first, InputIt last, Reduce f)->decltype(f(first, last)) {
    using OutputT = decltype(f(first, last));

    const size_t n = std::distance(first, last);
    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t batch = n / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    std::vector<OutputT>     fragments(thread_count);

    InputIt batch_last = first;

    size_t i = 0;
    for (; i != thread_count-1; ++i) {
        std::advance(batch_last, batch);
        threads[i] = std::thread(impl::reduce<InputIt, Reduce, OutputT>, 
            first, batch_last, f, std::ref(fragments[i]));
        std::advance(first, batch);
        batch_last = first;
    }
    impl::reduce(first, last, f, fragments[i]);

    for (auto &th : threads) th.join();
    return (fragments.size() == 1)
        ? fragments.front()
        : std::accumulate(fragments.begin()+1, fragments.end(), fragments.front());
}

}; //namespace bio

#endif