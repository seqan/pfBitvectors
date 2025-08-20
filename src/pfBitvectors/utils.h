//SPDX-FileCopyrightText: 2024 Simon Gene Gottlieb
//SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <array>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <vector>

namespace seqan::pfb {

template <size_t N>
inline std::array<std::bitset<N>, N+1> const leftshift_masks = []() {
    auto m = std::array<std::bitset<N>, N+1>{};
    m[0].flip();

    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] >> 1;
    }
    return m;
}();

template <size_t N>
inline std::array<std::bitset<N>, N+1> const rightshift_masks = []() {
    auto m = std::array<std::bitset<N>, N+1>{};
    m[0].flip();

    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] << 1;
    }
    return m;
}();

template <size_t N>
std::array<std::bitset<N>, N*2+1> const signed_rightshift_masks = []() {
    auto m = std::array<std::bitset<N>, N*2+1>{};
    m[0].flip();
    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] << 1;
    }
    for (size_t i{N+1}; i < N*2+1; ++i) {
        m[i] = m[i-1] << 1;
        m[i].set(0);
    }
    return m;
}();

template <size_t N>
size_t lshift_and_count(std::bitset<N> const& b, size_t shift) {
    auto const& mask = leftshift_masks<N>[shift];
    return (b & mask).count();
}

template <size_t N>
size_t rshift_and_count(std::bitset<N> const& b, size_t shift) {
    auto const& mask = rightshift_masks<N>[shift];
    return (b & mask).count();
}

template <size_t N>
size_t signed_rshift_and_count(std::bitset<N> const& b, size_t shift) {
    auto const& mask = signed_rightshift_masks<N>[shift];
    return (b & mask).count();
}

/**
 * an array of size N*2+1. This array masks the first x bits or last y bits, depending on the parameter
 *
 * 0: 1111
 * 1: 1110
 * 2: 1100
 * 3: 1000
 * 4: 0000
 * 5: 0001
 * 6: 0011
 * 7: 0111
 * 8: 1111
 */
template <size_t N>
alignas((N==512)?64:alignof(std::array<std::bitset<N>, N*2+1>))
std::array<std::bitset<N>, N*2+1> const skip_first_or_last_n_bits_masks = []() {
    auto m = std::array<std::bitset<N>, N*2+1>{};
    m[0].flip();
    for (size_t i{1}; i < N+1; ++i) {
        m[i] = m[i-1] << 1;
    }
    for (size_t i{N+1}; i < N*2+1; ++i) {
        m[i] = m[i-1] << 1;
        m[i].set(0);
    }
    return m;
}();

/*
 * e.g.: N=8
 * idx 0: count all 8 bits
 * idx 1: skip first bit and count tailing 7bits
 * ....
 * idx 7: skip 7 bits and count only last bit
 * idx 8: count zero bits
 * idx 9: count first bit and skip last 7bits
 * ...
 * idx 15: count first 7 bits and skip last bit
 * idx 16: count all bits
 */
template <size_t N>
size_t skip_first_or_last_n_bits_and_count(std::bitset<N> const& b, size_t idx) {
    auto const& mask = skip_first_or_last_n_bits_masks<N>[idx];
    return (b & mask).count();
}
}
