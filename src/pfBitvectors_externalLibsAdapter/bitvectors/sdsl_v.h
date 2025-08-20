// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <sdsl/bit_vectors.hpp>

namespace seqan::pfb {

struct SDSL_V {
    sdsl::bit_vector bitvector;
    sdsl::bit_vector::rank_1_type bv;

    SDSL_V() = default;
    SDSL_V(SDSL_V&&) = default;

    // constructor accepting view to bools
    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, bool>
    SDSL_V(range_t&& _range)
        : bitvector{[&]() {
            auto bv = sdsl::bit_vector(_range.size(), 0);
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
        , bv{&bitvector}
    {
    }

    auto symbol(size_t i) const -> bool {
        return bitvector[i];
    }

    auto rank(size_t i) const -> size_t {
        return bv.rank(i);
    }

    auto size() const -> size_t {
        return bv.size();
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, bv);
    }
};

}
