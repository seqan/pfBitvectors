<!--
    SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Paired-Blocks Flattened Bitvectors

A C++ Library that provides implementation of the "flattened bit vectors" and "Paired-Block bit vectors".


## Introduction
A common challenge is having (very long) strings and wanting to know the `rank` of a character `c` at a certain position.
The `rank` returns the numbers `c` in the first i position.

These string benefit from the knowing that only a very small alphabet is being used. (Alphabet e.g.: Only ACGT (DNA data) or only 26 characters).


## Available structures
We provide multiple implementation for bit vectors and strings with rank support.
Following classes provide bit vectors with rank support:
- `seqan::pfb::Bitvector<...>`
- `seqan::pfb::PairedBitvector<...>`

Following classes provide strings with rank support
- `seqan::pfb::FlattenedBitvectors2L<...>`
- `seqan::pfb::PairedFlattenedBitvectors2L<...>`


## Usage
### Setup
There are multiple ways of getting this up and going:
- via [CPM](https://github.com/cpm-cmake/CPM.cmake) (recommended): `CPMAddPackage("gh:seqan/pfBitvectors@1.0.0")`
- via [add_subdirectory](https://cmake.org/cmake/help/latest/command/add_subdirectory.html): `add_subdirectory(path/to/pfBitvectors EXCLUDE_FROM_ALL SYSTEM)`
- via gcc/clang: `-isystem /path/to/pfBitvectors/src`
- via headers only: copy the folder `src/pfBitvectors` to whereever you like

For any CMake related approach, you must link against `seqan::pbf` as `target_link_libraries(yourTarget .... seqan::pbf)`

### Example 1

An example to use a bit vector. Complete example can be seen in file [exampleBitvectors/main.cpp](examples/exampleBitvectors/main.cpp).
```c++
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
```

### Example 2
Another example which uses strings with rank support.  Complete example can be seen in file [exampleStrings/main.cpp](examples/exampleStrings/main.cpp).
```c++
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
```


## Citation
For academic work please cite:
- *Engineering rank queries on bit vectors and strings*; Simon Gene Gottlieb, Knut Reinert; preprint

## Authorship and Copyright
*pfBitvectors* is being developed by [Simon Gene Gottlieb](mailto.simon@gottliebtfreitag.de) at the [Algorithmic Bioinformatics Group](https://www.mi.fu-berlin.de/en/inf/groups/abi/index.html) of the [Freie University Berlin](https://www.fu-berlin.de/).
The *pfBitvectors* library is considered to be part of the [SeqAn project](https://www.seqan.de/).

## License
All files of this repository carry either a SPDX header or if not possible a .license file with the same information.
The library source code is under the BSD-3 license while code snippets or example are CC0-1.0.
The documentation (excluding its code snippets) has a CC-BY-4.0 license.
