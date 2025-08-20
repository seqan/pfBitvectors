// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "Bitvector1L.h"
#include "Bitvector2L.h"

namespace seqan::pfb {

template <size_t... Ns>
struct Bitvector {
    static_assert(sizeof...(Ns) != 0, "at least one level must be given");
    static_assert(sizeof...(Ns) <= 2, "to many levels, only one ore two level are possible");
};

template <size_t L0>
struct Bitvector<L0> : Bitvector1L<L0> {
    template <typename Archive>
    void serialize(Archive& ar) {
        Bitvector1L<L0>::serialize(ar);
    }
};

template <size_t L1, size_t L0>
struct Bitvector<L1, L0> : Bitvector2L<L1, L0> {
    template <typename Archive>
    void serialize(Archive& ar) {
        Bitvector2L<L1, L0>::serialize(ar);
    }
};


//template <size_t L1, size_t L0> : Bitvector2L<L1, L
}
