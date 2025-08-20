// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <bitset>
#include <cstddef>

#include "utils.h"

namespace seqan::pfb {

// helper to compute 'alignas' values
constexpr inline auto alignAsValue(size_t bits) -> size_t {
    if (bits % 512 == 0) return 64;
    if (bits % 256 == 0) return 32;
    if (bits % 128 == 0) return 16;
    if (bits %  64 == 0)  return 8;
    return 1;
}
// helper to compute 'alignas' values
constexpr inline auto minAlignAsValue(auto... bits) {
    return std::min(alignAsValue(bits)...);
}


template <size_t N, typename Archive>
void loadBV(std::bitset<N>& b, Archive& ar) {
    b = std::bitset<N>{};
    for (size_t i{0}; i < N; i += 64) {
        uint64_t v{};
        ar(v);
        b = b | (std::bitset<N>{v} << i);
    }
}
template <size_t N, typename Archive>
void saveBV(std::bitset<N> const& b, Archive& ar) {
    (void)b;
    (void)ar;
    static constexpr auto mask = std::bitset<N>{~uint64_t{0}};

    // saving in 64bit blocks
    for (size_t i{0}; i < N; i += 64) {
        auto v = ((b >> i) & mask).to_ullong();
        ar(v);
    }
}
template <size_t N, bool Align=true>
struct alignas(std::max(alignof(std::bitset<N>), Align?alignAsValue(N):size_t{1})) AlignedBitset {
    std::bitset<N> bits;

    decltype(auto) operator[](size_t i) {
        return bits[i];
    }
    decltype(auto) operator[](size_t i) const {
        return bits[i];
    }

    auto count() const {
        return bits.count();
    }
    auto size() const -> size_t {
        return N;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        saveBV(bits, ar);
    }

    template <typename Archive>
    void load(Archive& ar) {
        loadBV(bits, ar);
    }
};

}
