<!--
    SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Paired-Blocks Flattened Bitvectors

A C++ Library that provides implementation of the "flattened bit vectors" and "Paired-Block bit vectors".


## Introduction
A common challenge is having (very long) strings and requesting wanting to know the `rank` of a character `c` at a certain position.
The `rank` returns the numbers `c` in the first i position.

These string benefit from the knowing that only a very small alphabet is being used. (Alphabet e.g.: Only ACGT (DNA data) or only 26 characters).


## Available structures
We provide 3 structures/classes:
- `PairedBitvector`: bit vector with rank support
- `FlattenedBitvectors`: string with rank support
- `PairedFlattenedBitvectors`: string with rank support

## Usage
!TODO write the rest


## Citation
For academic work please cite:
- *Engineering rank queries on bit vectors and strings*; Simon Gene Gottlieb, Knut Reinert; preprint

## Authorship and Copyright
*pfBitvectors* is being developed by [Simon Gene Gottlieb](mailto.simon@gottliebtfreitag.de) at the [Algorithmic Bioinformatics Group](https://www.mi.fu-berlin.de/en/inf/groups/abi/index.html) of the [Freie University Berlin](https://www.fu-berlin.de/).
The *pfBitvectors* library is considered to be part of the [SeqAn project](https://www.seqan.de/).

## License
All files of this repository carries either have SPDX header or if not possible a .license file lies next to them.
The library source code is under the BSD-3 license while code snippets or example are CC0-1.0.
The documentation has a CC-BY-4.0 license.
