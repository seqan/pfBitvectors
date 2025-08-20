// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <rank9.h>

namespace seqan::pfb {

struct Rank9 {
    std::vector<uint64_t> bitvector;
    rank9 bv;
    size_t totalSize;

    Rank9() = default;
    Rank9(Rank9&&) = default;

    // constructor accepting view to bools
    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, bool>
    Rank9(range_t&& _range)
        : bitvector{[&_range]() {
            auto bitvector = std::vector<uint64_t>(_range.size(), 0);
            bitvector.resize(_range.size()/64 + 1);
            auto iter = _range.begin();
            size_t i = 0;
            for (;i < _range.size() && iter != _range.end(); ++i, ++iter) {
                auto id = i / 64;
                auto offset = i % 64;
                bitvector[id] |= (static_cast<uint64_t>(*iter)<<offset);
            }
            if (i != _range.size() || iter != _range.end()) {
                throw std::runtime_error{"Error constructing FlatRank (should not happen, please report to seqan::pfb library."};
            }
            return bitvector;
        }()}
        , bv{bitvector.data(), _range.size()+1}
        , totalSize{_range.size()}

    {
    }

    auto symbol(size_t i) const -> bool {
        auto id = i / 64;
        auto offset = i % 64;
        return (bitvector[id] >> offset) & 1;
    }

    auto rank(size_t i) const -> size_t {
        //hack, since rank9 is not const correct
        return const_cast<rank9&>(bv).rank(i);
    }

    auto size() const -> size_t {
        return totalSize;
    }

    auto space_usage() const -> size_t {
        return 0;
//        bitvector.size() * 8 + rank9.
    }
};

}
