// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <cereal/archives/binary.hpp>
#include <cstdlib>
#include <pfBitvectors/pfBitvectors.h>
#include <pfBitvectors_externalLibsAdapter/all.h>
#include <pfBitvectors_test_utils/utils.h>
#include <sstream>

using AllStrings = Variant<
    Instance<seqan::pfb::FlattenedBitvectors2L,   64,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  128,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  256,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  512,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L, 1024,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L, 2048,  4096>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,   64, 65536>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  128, 65536>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  256, 65536>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L,  512, 65536>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L, 1024, 65536>::Type,
    Instance<seqan::pfb::FlattenedBitvectors2L, 2048, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L,   64, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L,  128, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L,  256, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L,  512, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L, 1024, 65536>::Type,
    Instance<seqan::pfb::PairedFlattenedBitvectors2L, 2048, 65536>::Type,
    seqan::pfb::MultiBitvector,
#ifdef PFBITVECTORS_USE_SDSL
    seqan::pfb::Sdsl_wt_bldc,
    seqan::pfb::Sdsl_wt_epr,
#endif
    Delimiter /*delimiter, is ignored*/
>;


template <size_t min, size_t range>
auto generateText(size_t length) -> std::vector<uint8_t> {
    auto seed = []() -> size_t {
        if (auto ptr = std::getenv("FPB_SEED")) {
            return std::stoull(ptr);
        }
        return 0;
    }();
    srand(seed);

    auto text = std::vector<uint8_t>{};
    for (size_t i{0}; i<length; ++i) {
        text.push_back((rand() % range) + min);
    }
    return text;
}

TEST_CASE("check if rank on the symbol vectors is working, all sizes", "[string][all_sizes]") {
    auto testSigma = []<size_t Sigma>() {
        INFO("Sigma " << Sigma);
        call_with_templates([&]<template <size_t> typename _String>() {
            using String = _String<Sigma>;
            auto vector_name = getName<String>();
            INFO(vector_name);

            auto text = generateText<0, String::Sigma>(1000);

            auto vec = String{std::span{text}};
            REQUIRE(vec.size() == text.size());
            {
                for (size_t i{0}; i < text.size(); ++i) {
                    INFO(i);
                    CHECK(vec.symbol(i) == text.at(i));
                }
            }
            {
                auto countRank = [&](size_t idx, uint8_t sym) {
                    size_t acc{};
                    for (size_t i{0}; i < idx; ++i) {
                        acc = acc + (text[i] == sym);
                    }
                    return acc;
                };
                for (size_t i{0}; i <= text.size(); ++i) {
                    INFO(i);
                    for (size_t symb{}; symb < String::Sigma; ++symb) {
                        CHECK(vec.rank(i, symb) == countRank(i, symb));
                    }
                }
            }
        }, AllStrings{});
    };

    SECTION("test different sizes of alphabets") {
        testSigma.operator()<4>();
        testSigma.operator()<5>();
        testSigma.operator()<6>();
        testSigma.operator()<21>();
        testSigma.operator()<255>();
    }
}

TEST_CASE("hand counted, test with 255 alphabet", "[string][255][small]") {

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};

    SECTION("checks") {
        call_with_templates([&]<template <size_t> typename _String>() {
            using String = _String<255>;
            auto vector_name = getName<String>();
            INFO(vector_name);

            auto vec = String{std::span{text}};
            // checking construction size
            REQUIRE(vec.size() == text.size());

            // checking symbols
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vec.symbol(i) == text.at(i));
            }
            // test complete vector on symbol()"
            CHECK(vec.symbol( 0) == 'H');
            CHECK(vec.symbol( 1) == 'a');
            CHECK(vec.symbol( 2) == 'l');
            CHECK(vec.symbol( 3) == 'l');
            CHECK(vec.symbol( 4) == 'o');
            CHECK(vec.symbol( 5) == ' ');
            CHECK(vec.symbol( 6) == 'W');
            CHECK(vec.symbol( 7) == 'e');
            CHECK(vec.symbol( 8) == 'l');
            CHECK(vec.symbol( 9) == 't');

            // test complete vector on rank()
            CHECK(vec.rank( 0, ' ') == 0);
            CHECK(vec.rank( 1, ' ') == 0);
            CHECK(vec.rank( 2, ' ') == 0);
            CHECK(vec.rank( 3, ' ') == 0);
            CHECK(vec.rank( 4, ' ') == 0);
            CHECK(vec.rank( 5, ' ') == 0);
            CHECK(vec.rank( 6, ' ') == 1);
            CHECK(vec.rank( 7, ' ') == 1);
            CHECK(vec.rank( 8, ' ') == 1);
            CHECK(vec.rank( 9, ' ') == 1);
            CHECK(vec.rank(10, ' ') == 1);

            CHECK(vec.rank( 0, 'H') == 0);
            CHECK(vec.rank( 1, 'H') == 1);
            CHECK(vec.rank( 2, 'H') == 1);
            CHECK(vec.rank( 3, 'H') == 1);
            CHECK(vec.rank( 4, 'H') == 1);
            CHECK(vec.rank( 5, 'H') == 1);
            CHECK(vec.rank( 6, 'H') == 1);
            CHECK(vec.rank( 7, 'H') == 1);
            CHECK(vec.rank( 8, 'H') == 1);
            CHECK(vec.rank( 9, 'H') == 1);
            CHECK(vec.rank(10, 'H') == 1);

            CHECK(vec.rank( 0, 'W') == 0);
            CHECK(vec.rank( 1, 'W') == 0);
            CHECK(vec.rank( 2, 'W') == 0);
            CHECK(vec.rank( 3, 'W') == 0);
            CHECK(vec.rank( 4, 'W') == 0);
            CHECK(vec.rank( 5, 'W') == 0);
            CHECK(vec.rank( 6, 'W') == 0);
            CHECK(vec.rank( 7, 'W') == 1);
            CHECK(vec.rank( 8, 'W') == 1);
            CHECK(vec.rank( 9, 'W') == 1);
            CHECK(vec.rank(10, 'W') == 1);

            CHECK(vec.rank( 0, 'a') == 0);
            CHECK(vec.rank( 1, 'a') == 0);
            CHECK(vec.rank( 2, 'a') == 1);
            CHECK(vec.rank( 3, 'a') == 1);
            CHECK(vec.rank( 4, 'a') == 1);
            CHECK(vec.rank( 5, 'a') == 1);
            CHECK(vec.rank( 6, 'a') == 1);
            CHECK(vec.rank( 7, 'a') == 1);
            CHECK(vec.rank( 8, 'a') == 1);
            CHECK(vec.rank( 9, 'a') == 1);
            CHECK(vec.rank(10, 'a') == 1);

            CHECK(vec.rank( 0, 'e') == 0);
            CHECK(vec.rank( 1, 'e') == 0);
            CHECK(vec.rank( 2, 'e') == 0);
            CHECK(vec.rank( 3, 'e') == 0);
            CHECK(vec.rank( 4, 'e') == 0);
            CHECK(vec.rank( 5, 'e') == 0);
            CHECK(vec.rank( 6, 'e') == 0);
            CHECK(vec.rank( 7, 'e') == 0);
            CHECK(vec.rank( 8, 'e') == 1);
            CHECK(vec.rank( 9, 'e') == 1);
            CHECK(vec.rank(10, 'e') == 1);

            CHECK(vec.rank( 0, 'l') == 0);
            CHECK(vec.rank( 1, 'l') == 0);
            CHECK(vec.rank( 2, 'l') == 0);
            CHECK(vec.rank( 3, 'l') == 1);
            CHECK(vec.rank( 4, 'l') == 2);
            CHECK(vec.rank( 5, 'l') == 2);
            CHECK(vec.rank( 6, 'l') == 2);
            CHECK(vec.rank( 7, 'l') == 2);
            CHECK(vec.rank( 8, 'l') == 2);
            CHECK(vec.rank( 9, 'l') == 3);
            CHECK(vec.rank(10, 'l') == 3);

            CHECK(vec.rank( 0, 'o') == 0);
            CHECK(vec.rank( 1, 'o') == 0);
            CHECK(vec.rank( 2, 'o') == 0);
            CHECK(vec.rank( 3, 'o') == 0);
            CHECK(vec.rank( 4, 'o') == 0);
            CHECK(vec.rank( 5, 'o') == 1);
            CHECK(vec.rank( 6, 'o') == 1);
            CHECK(vec.rank( 7, 'o') == 1);
            CHECK(vec.rank( 8, 'o') == 1);
            CHECK(vec.rank( 9, 'o') == 1);
            CHECK(vec.rank(10, 'o') == 1);

            CHECK(vec.rank( 0, 't') == 0);
            CHECK(vec.rank( 1, 't') == 0);
            CHECK(vec.rank( 2, 't') == 0);
            CHECK(vec.rank( 3, 't') == 0);
            CHECK(vec.rank( 4, 't') == 0);
            CHECK(vec.rank( 5, 't') == 0);
            CHECK(vec.rank( 6, 't') == 0);
            CHECK(vec.rank( 7, 't') == 0);
            CHECK(vec.rank( 8, 't') == 0);
            CHECK(vec.rank( 9, 't') == 0);
            CHECK(vec.rank(10, 't') == 1);

            // check all other characters have rank 0
            auto ignore = std::unordered_set<size_t>{' ', 'H', 'W', 'a', 'e', 'l', 'o', 't'};
            for (size_t s{0}; s < String::Sigma; ++s) {
                if (ignore.contains(s)) continue;
                for (size_t i{0}; i < 11; ++i) {
                    CHECK(vec.rank(i, s) == 0);
                }
            }

            // test complete vec 'H' for prefix_rank()
            CHECK(vec.prefix_rank( 0, ' ') == 0);
            CHECK(vec.prefix_rank( 1, ' ') == 0);
            CHECK(vec.prefix_rank( 2, ' ') == 0);
            CHECK(vec.prefix_rank( 3, ' ') == 0);
            CHECK(vec.prefix_rank( 4, ' ') == 0);
            CHECK(vec.prefix_rank( 5, ' ') == 0);
            CHECK(vec.prefix_rank( 6, ' ') == 0);
            CHECK(vec.prefix_rank( 7, ' ') == 0);
            CHECK(vec.prefix_rank( 8, ' ') == 0);
            CHECK(vec.prefix_rank( 9, ' ') == 0);
            CHECK(vec.prefix_rank(10, ' ') == 0);

            CHECK(vec.prefix_rank( 0, 'H') == 0);
            CHECK(vec.prefix_rank( 1, 'H') == 0);
            CHECK(vec.prefix_rank( 2, 'H') == 0);
            CHECK(vec.prefix_rank( 3, 'H') == 0);
            CHECK(vec.prefix_rank( 4, 'H') == 0);
            CHECK(vec.prefix_rank( 5, 'H') == 0);
            CHECK(vec.prefix_rank( 6, 'H') == 1);
            CHECK(vec.prefix_rank( 7, 'H') == 1);
            CHECK(vec.prefix_rank( 8, 'H') == 1);
            CHECK(vec.prefix_rank( 9, 'H') == 1);
            CHECK(vec.prefix_rank(10, 'H') == 1);

            CHECK(vec.prefix_rank( 0, 'W') == 0);
            CHECK(vec.prefix_rank( 1, 'W') == 1);
            CHECK(vec.prefix_rank( 2, 'W') == 1);
            CHECK(vec.prefix_rank( 3, 'W') == 1);
            CHECK(vec.prefix_rank( 4, 'W') == 1);
            CHECK(vec.prefix_rank( 5, 'W') == 1);
            CHECK(vec.prefix_rank( 6, 'W') == 2);
            CHECK(vec.prefix_rank( 7, 'W') == 2);
            CHECK(vec.prefix_rank( 8, 'W') == 2);
            CHECK(vec.prefix_rank( 9, 'W') == 2);
            CHECK(vec.prefix_rank(10, 'W') == 2);

            CHECK(vec.prefix_rank( 0, 'a') == 0);
            CHECK(vec.prefix_rank( 1, 'a') == 1);
            CHECK(vec.prefix_rank( 2, 'a') == 1);
            CHECK(vec.prefix_rank( 3, 'a') == 1);
            CHECK(vec.prefix_rank( 4, 'a') == 1);
            CHECK(vec.prefix_rank( 5, 'a') == 1);
            CHECK(vec.prefix_rank( 6, 'a') == 2);
            CHECK(vec.prefix_rank( 7, 'a') == 3);
            CHECK(vec.prefix_rank( 8, 'a') == 3);
            CHECK(vec.prefix_rank( 9, 'a') == 3);
            CHECK(vec.prefix_rank(10, 'a') == 3);

            CHECK(vec.prefix_rank( 0, 'e') == 0);
            CHECK(vec.prefix_rank( 1, 'e') == 1);
            CHECK(vec.prefix_rank( 2, 'e') == 2);
            CHECK(vec.prefix_rank( 3, 'e') == 2);
            CHECK(vec.prefix_rank( 4, 'e') == 2);
            CHECK(vec.prefix_rank( 5, 'e') == 2);
            CHECK(vec.prefix_rank( 6, 'e') == 3);
            CHECK(vec.prefix_rank( 7, 'e') == 4);
            CHECK(vec.prefix_rank( 8, 'e') == 4);
            CHECK(vec.prefix_rank( 9, 'e') == 4);
            CHECK(vec.prefix_rank(10, 'e') == 4);

            CHECK(vec.prefix_rank( 0, 'l') == 0);
            CHECK(vec.prefix_rank( 1, 'l') == 1);
            CHECK(vec.prefix_rank( 2, 'l') == 2);
            CHECK(vec.prefix_rank( 3, 'l') == 2);
            CHECK(vec.prefix_rank( 4, 'l') == 2);
            CHECK(vec.prefix_rank( 5, 'l') == 2);
            CHECK(vec.prefix_rank( 6, 'l') == 3);
            CHECK(vec.prefix_rank( 7, 'l') == 4);
            CHECK(vec.prefix_rank( 8, 'l') == 5);
            CHECK(vec.prefix_rank( 9, 'l') == 5);
            CHECK(vec.prefix_rank(10, 'l') == 5);

            CHECK(vec.prefix_rank( 0, 'o') == 0);
            CHECK(vec.prefix_rank( 1, 'o') == 1);
            CHECK(vec.prefix_rank( 2, 'o') == 2);
            CHECK(vec.prefix_rank( 3, 'o') == 3);
            CHECK(vec.prefix_rank( 4, 'o') == 4);
            CHECK(vec.prefix_rank( 5, 'o') == 4);
            CHECK(vec.prefix_rank( 6, 'o') == 5);
            CHECK(vec.prefix_rank( 7, 'o') == 6);
            CHECK(vec.prefix_rank( 8, 'o') == 7);
            CHECK(vec.prefix_rank( 9, 'o') == 8);
            CHECK(vec.prefix_rank(10, 'o') == 8);

            CHECK(vec.prefix_rank( 0, 't') == 0);
            CHECK(vec.prefix_rank( 1, 't') == 1);
            CHECK(vec.prefix_rank( 2, 't') == 2);
            CHECK(vec.prefix_rank( 3, 't') == 3);
            CHECK(vec.prefix_rank( 4, 't') == 4);
            CHECK(vec.prefix_rank( 5, 't') == 5);
            CHECK(vec.prefix_rank( 6, 't') == 6);
            CHECK(vec.prefix_rank( 7, 't') == 7);
            CHECK(vec.prefix_rank( 8, 't') == 8);
            CHECK(vec.prefix_rank( 9, 't') == 9);
            CHECK(vec.prefix_rank(10, 't') == 9);

            CHECK(vec.prefix_rank( 0, 'z') ==  0);
            CHECK(vec.prefix_rank( 1, 'z') ==  1);
            CHECK(vec.prefix_rank( 2, 'z') ==  2);
            CHECK(vec.prefix_rank( 3, 'z') ==  3);
            CHECK(vec.prefix_rank( 4, 'z') ==  4);
            CHECK(vec.prefix_rank( 5, 'z') ==  5);
            CHECK(vec.prefix_rank( 6, 'z') ==  6);
            CHECK(vec.prefix_rank( 7, 'z') ==  7);
            CHECK(vec.prefix_rank( 8, 'z') ==  8);
            CHECK(vec.prefix_rank( 9, 'z') ==  9);
            CHECK(vec.prefix_rank(10, 'z') == 10);

            // check all_ranks() is equal to prefix_rank() and rank()
            for (size_t idx{0}; idx < vec.size(); ++idx) {
                auto [rank, prefix] = vec.all_ranks_and_prefix_ranks(idx);
                auto rank2 = vec.all_ranks(idx);
                for (size_t symb{1}; symb < String::Sigma; ++symb) {
                    INFO(idx);
                    INFO(symb);
                    CHECK(rank[symb] == vec.rank(idx, symb));
                    CHECK(rank2[symb] == vec.rank(idx, symb));
                    CHECK(prefix[symb] == vec.prefix_rank(idx, symb));
                }
            }
        }, AllStrings{});
    }
}

TEST_CASE("check symbol vectors construction on text longer than 255 characters", "[string][255][large]") {
    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                     254, 254, 254, 254, 254, 254, 254, 254, 254, 254,
                                    };

    SECTION("checks") {
        call_with_templates([&]<template <size_t> typename _String>() {
            using String = _String<255>;

            auto vector = String{text};

            REQUIRE(vector.size() == text.size());

            // check that symbol() call works
            for (size_t i{0}; i < text.size(); ++i) {
                INFO(i);
                CHECK(vector.symbol(i) == text.at(i));
            }

            auto countRank = [&](size_t idx, uint8_t sym) {
                size_t acc{};
                for (size_t i{0}; i < idx; ++i) {
                    acc = acc + (text[i] == sym);
                }
                return acc;
            };
            auto countPrefixRank = [&](size_t idx, uint8_t sym) {
                size_t acc{};
                for (size_t i{0}; i < idx; ++i) {
                    acc = acc + (text[i] < sym);
                }
                return acc;
            };


            // check all_ranks() is equal to prefix_rank() and rank()
            for (size_t idx{0}; idx <= vector.size(); ++idx) {
                auto [rank, prefix] = vector.all_ranks_and_prefix_ranks(idx);
                auto rank2 = vector.all_ranks(idx);
                for (size_t symb{1}; symb < String::Sigma; ++symb) {
                    INFO(idx);
                    INFO(symb);
                    CHECK(countRank(idx, symb) == vector.rank(idx, symb));
                    CHECK(countPrefixRank(idx, symb) == vector.prefix_rank(idx, symb));
                    vector.rank(idx, symb);
                    vector.prefix_rank(idx, symb);
                    CHECK(rank[symb] == vector.rank(idx, symb));
                    CHECK(rank2[symb] == vector.rank(idx, symb));
                    CHECK(countRank(idx, symb) == rank[symb]);
                    CHECK(countPrefixRank(idx, symb) == prefix[symb]);
                    CHECK(prefix[symb] == vector.prefix_rank(idx, symb));
                }
            }
        }, AllStrings{});
    }
}
