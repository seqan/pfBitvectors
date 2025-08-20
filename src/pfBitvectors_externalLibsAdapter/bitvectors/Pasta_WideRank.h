// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/rank.hpp>
#include <pasta/bit_vector/support/wide_rank.hpp>

#include <ranges>

namespace seqan::pfb {

struct WideRank {
    pasta::BitVector bv;
    pasta::WideRank<pasta::OptimizedFor::DONT_CARE, pasta::BitVector> rs{bv};

    WideRank() = default;
    WideRank(WideRank&&) = default;

    // constructor accepting view to bools
    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, bool>
    WideRank(range_t&& _range)
        : bv{[&]() {
            pasta::BitVector bv{_range.size(), 0};
            auto iter = _range.begin();
            size_t i = 0;
            for (;i < bv.size() && iter != _range.end(); ++i, ++iter) {
                bv[i] = *iter;
            }
            if (i != bv.size() || iter != _range.end()) {
                throw std::runtime_error{"Error constructing FlatRank (should not happen, please report to seqan::pfb library."};
            }
            return bv;
        }()}
    {
    }

    auto symbol(size_t i) const -> bool {
        return bv[i];
    }

    auto rank(size_t i) const -> size_t {
        return rs.rank1(i);
    }

    auto size() const -> size_t {
        return bv.size();
    }

    auto space_usage() const -> size_t {
        return bv.space_usage() + rs.space_usage();
    }
};

}
