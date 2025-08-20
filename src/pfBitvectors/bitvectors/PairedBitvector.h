// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "PairedBitvector1L.h"
#include "PairedBitvector2L.h"

namespace seqan::pfb {

template <size_t... Ns>
struct PairedBitvector {
    static_assert(sizeof...(Ns) != 0, "at least one level must be given");
    static_assert(sizeof...(Ns) <= 2, "to many levels, only one ore two level are possible");
};

template <size_t L0>
struct PairedBitvector<L0> : PairedBitvector1L<L0> {
    template <typename Archive>
    void serialize(Archive& ar) {
        PairedBitvector1L<L0>::serialize(ar);
    }
};

template <size_t L1, size_t L0>
struct PairedBitvector<L1, L0> : PairedBitvector2L<L1, L0> {
    template <typename Archive>
    void serialize(Archive& ar) {
        PairedBitvector2L<L1, L0>::serialize(ar);
    }
};


//template <size_t L1, size_t L0> : Bitvector2L<L1, L
}
