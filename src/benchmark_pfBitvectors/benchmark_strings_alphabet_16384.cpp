// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "BenchSize.h"

#include <catch2/catch_all.hpp>
#include <cereal/archives/binary.hpp>
#include <cstddef>
#include <nanobench.h>
#include <pfBitvectors/pfBitvectors.h>
#include <pfBitvectors_externalLibsAdapter/all.h>
#include <pfBitvectors_test_utils/utils.h>
#include <string>

namespace {
    #define SIGMA 16384
    #define STRINGIFY(x) #x
    #define TOSTRING(x) STRINGIFY(x)
    constexpr static size_t Sigma = SIGMA;
    #define SIGMA_STR TOSTRING(SIGMA)
}

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
    Delimiter /*delimiter, is ignored*/
>;


template <size_t min, size_t range>
auto generateText(size_t length) -> std::vector<uint8_t> {
    auto rng = ankerl::nanobench::Rng{};

    auto text = std::vector<uint8_t>{};
    for (size_t i{0}; i<length; ++i) {
        text.push_back(rng.bounded(range) + min);
    }
    return text;
}

template <size_t min, size_t range>
auto generateText() -> std::vector<uint8_t> const& {
    static auto text = []() -> std::vector<uint8_t> {
        auto rng = ankerl::nanobench::Rng{};

        // generates string with values between 1-4
        auto size = []() -> size_t {
            auto ptr = std::getenv("STRINGSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 1'000'000;
            #else
                return 1'000;
            #endif
        }();
        return generateText<min, range>(size);
    }();
    return text;
}


TEST_CASE("benchmark strings c'tor operation - " SIGMA_STR " alphabet", "[string][" SIGMA_STR "][time][ctor]") {
    auto const& text = generateText<0, Sigma>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            bench.run(name, [&]() {
                auto str = String{text};
                ankerl::nanobench::doNotOptimizeAway(const_cast<String const&>(str));
            });
        }, AllStrings{});
    }
}

TEST_CASE("benchmark vectors symbol() operations - " SIGMA_STR " alphabet", "[string][" SIGMA_STR "][time][symbol]") {
    auto const& text = generateText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("symbol()")
             .relative(true)
             .batch(text.size());

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllStrings{});
    }
}

TEST_CASE("benchmark vectors rank() operations - " SIGMA_STR " alphabet", "[string][" SIGMA_STR "][time][rank]") {
    auto const& text = generateText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("rank()")
             .relative(true);

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.rank(rng.bounded(text.size()+1), rng.bounded(Sigma));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllStrings{});
    }
}

TEST_CASE("benchmark vectors prefix_rank() operations - " SIGMA_STR " alphabet", "[string][" SIGMA_STR "][time][prefix_rank]") {
    auto const& text = generateText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("prefix_rank()")
             .relative(true);

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.prefix_rank(rng.bounded(text.size()+1), rng.bounded(Sigma));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllStrings{});
    }
}

TEST_CASE("benchmark vectors all_ranks() operations - " SIGMA_STR " alphabet", "[string][" SIGMA_STR "][time][all_ranks]") {
    auto const& text = generateText<0, Sigma>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks()")
             .relative(true);

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto rng = ankerl::nanobench::Rng{};

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.all_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllStrings{});
    }
}

TEST_CASE("benchmark vectors all_ranks_and_prefix_ranks() operations - Sigma alphabet", "[string][" SIGMA_STR "][time][all_ranks_and_prefix_ranks]") {
    auto const& text = generateText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("all_ranks_and_prefix_ranks()")
             .relative(true);

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};

            bench.run(name, [&]() {
                auto v = str.all_ranks_and_prefix_ranks(rng.bounded(text.size()+1));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllStrings{});
    }
}
TEST_CASE("benchmark vectors in size - alphabet " SIGMA_STR, "[string][" SIGMA_STR "][size]") {
    auto const& text = generateText<0, Sigma>();
    auto rng = ankerl::nanobench::Rng{};

    SECTION("benchmarking") {
        BenchSize benchSize;
        benchSize.baseSize = std::ceil(std::log2(Sigma));
        benchSize.entries[0][2] = "bits/char";
        benchSize.entries[0][4] = "alphabet " SIGMA_STR;

        call_with_templates([&]<template <size_t> class _String>() {
            using String = _String<Sigma>;
            auto name = getName<String>();
            INFO(name);

            auto str = String{text};
            auto size = [&]() {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(str);
                return ofs.str().size();
            }();
            benchSize.addEntry({
                .name = name,
                .size = size,
                .text_size = text.size(),
                .bits_per_char = (size*8)/double(text.size())
            });
        }, AllStrings{});
    }
}
