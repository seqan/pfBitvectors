// SPDX-FileCopyrightText: 2025 Gottlieb+Freitag <info@gottliebtfreitag.de>
// SPDX-License-Identifier: CC0-1.0

#include <iostream>
#include <pfBitvectors/pfBitvectors.h>
#include <vector>
int main() {
    std::vector<uint8_t> values{0, 1, 2, 1, 0, 1, 2, 1, 2};
    // String with rank support, is using 3 different character values, blocks/superblocks are of size 512 and 65536
    auto string = seqan::pfb::FlattenedBitvectors2L<3, 512, 65536>{values};

    std::cout << "the first 4 chars contain " << string.rank(4, 0) << " characters of value 0\n";
    std::cout << "the character at index 5 has value " << string.symbol(5) << "\n";
    std::cout << "the string has length " << string.size() << "\n";
}
