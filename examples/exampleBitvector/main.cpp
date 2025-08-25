// SPDX-FileCopyrightText: 2025 Gottlieb+Freitag <info@gottliebtfreitag.de>
// SPDX-License-Identifier: CC0-1.0

#include <iostream>
#include <pfBitvectors/pfBitvectors.h>
#include <vector>
int main() {
    std::vector<bool> values{true, false, true, false};
    // A bit vector that uses 512-bit on the lowest level, and 65536-bits on the next one
    auto bitvector = seqan::pfb::Bitvector<512, 65536>{values};

    std::cout << "the first 3 bits have " << bitvector.rank(3) << " ones\n";
    std::cout << "the bit with index 3 has the value " << bitvector.symbol(3)\n";
    std::cout << "the bitvector is of length " << bitvector.size() <<"\n";
}
