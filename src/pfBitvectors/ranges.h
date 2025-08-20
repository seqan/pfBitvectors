// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <ranges>
#include <seqan-std/chunk_view.hpp>

namespace seqan::pfb {

constexpr inline auto view_bool_as_uint64 = seqan::stl::views::chunk(64) | std::views::transform([](auto r) -> uint64_t {
    auto v = uint64_t{};
    auto iter = r.begin();
    for (size_t j{0}; iter != r.end(); ++j, ++iter) {
        v = v | (uint64_t{*iter} << j);
    }
    return v;
});

template <size_t N>
constexpr inline auto view_as_bitset = seqan::stl::views::chunk(N / 64) | std::views::transform([](auto r) -> std::bitset<N> {
    static_assert(N % 64 == 0, "must be a multiple of 64");
    auto v = std::bitset<N>{};
    auto iter = r.begin();
    for (size_t j{0}; iter != r.end(); ++j, ++iter) {
        v = v | (std::bitset<N>{*iter} << (64*j));
    }
    return v;
});

}
