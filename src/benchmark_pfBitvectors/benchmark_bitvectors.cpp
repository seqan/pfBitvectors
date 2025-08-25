// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "BenchSize.h"

#include <catch2/catch_all.hpp>
#include <cereal/archives/binary.hpp>
#include <cstdlib>
#include <pfBitvectors/pfBitvectors.h>
#include <pfBitvectors_externalLibsAdapter/all.h>
#include <pfBitvectors_test_utils/utils.h>
#include <fstream>
#include <nanobench.h>

using AllBitvectors = std::variant<
#ifdef PFBITVECTORS_USE_PASTA
    seqan::pfb::FlatRank,
    seqan::pfb::WideRank,
#endif
#ifdef PFBITVECTORS_USE_SDSL
    seqan::pfb::SDSL_V,
    seqan::pfb::SDSL_V5,
#endif
#ifdef PFBITVECTORS_USE_SUX
    seqan::pfb::Rank9,
#endif
// Single layer bitvectors
//    seqan::pfb::Bitvector<  64>,
//    seqan::pfb::Bitvector< 128>,
//    seqan::pfb::Bitvector< 256>,
//    seqan::pfb::Bitvector< 512>,
//    seqan::pfb::Bitvector<1024>,
//    seqan::pfb::Bitvector<2048>,
//    seqan::pfb::PairedBitvector<  64>,
//    seqan::pfb::PairedBitvector< 128>,
//    seqan::pfb::PairedBitvector< 256>,
//    seqan::pfb::PairedBitvector< 512>,
//    seqan::pfb::PairedBitvector<1024>,
//    seqan::pfb::PairedBitvector<2048>,

// Two layer bitvectors
    seqan::pfb::Bitvector<  64, 65536>,
//    seqan::pfb::Bitvector< 128, 65536>,
//    seqan::pfb::Bitvector< 256, 65536>,
    seqan::pfb::Bitvector< 512, 65536>,
//    seqan::pfb::Bitvector<1024, 65536>,
//    seqan::pfb::Bitvector<2048, 65536>,
    seqan::pfb::PairedBitvector<  64, 65536>,
//    seqan::pfb::PairedBitvector< 128, 65536>,
//    seqan::pfb::PairedBitvector< 256, 65536>,
    seqan::pfb::PairedBitvector< 512, 65536>,
//    seqan::pfb::PairedBitvector<1024, 65536>,
//    seqan::pfb::PairedBitvector<2048, 65536>,
    std::monostate /*delimiter, is ignored*/
>;

namespace {
auto generateText() -> std::vector<bool> const& {
    static auto text = []() -> std::vector<bool> {
        auto rng = ankerl::nanobench::Rng{};

        auto size = []() -> size_t {
            auto ptr = std::getenv("BITVECTORSIZE");
            if (ptr) {
                return std::stoull(ptr);
            }
            #ifdef NDEBUG
                return 10'000'000;
            #else
                return 100'000;
            #endif
        }();

        auto text = std::vector<bool>{};
        for (size_t i{0}; i<size; ++i) {
            text.push_back(rng.bounded(4) == 0);
        }
        return text;
    }();
    return text;
}
}

TEST_CASE("benchmark bit vectors ctor run times", "[bitvector][time][ctor]") {
    auto bench_ctor = ankerl::nanobench::Bench{};
    bench_ctor.title("c'tor()")
              .relative(true);

    auto& text = generateText();

    SECTION("benchmarking") {
        call_with_templates([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            bench_ctor.batch(text.size()).run(vector_name, [&]() {
                auto vec = Vector{text};
                ankerl::nanobench::doNotOptimizeAway(vec.rank(0));
            });
        }, AllBitvectors{});
    }
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[bitvector][time][symbol]") {

    auto& text = generateText();

    SECTION("benchmarking - symbol") {
        auto bench_symbol = ankerl::nanobench::Bench{};
        bench_symbol.title("symbol()")
                    .relative(true);

        bench_symbol.epochs(10);
        bench_symbol.minEpochTime(std::chrono::milliseconds{10});
        call_with_templates([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench_symbol.run(vector_name, [&]() {
                auto v = vec.symbol(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllBitvectors{});
    }
}

TEST_CASE("benchmark bit vectors rank and symbol run times", "[bitvector][time][rank]") {

    auto& text = generateText();

    SECTION("benchmarking - rank") {
        auto bench_rank = ankerl::nanobench::Bench{};
        bench_rank.title("rank()")
                  .relative(true);

        bench_rank.epochs(20);
        bench_rank.minEpochTime(std::chrono::milliseconds{10});
        bench_rank.minEpochIterations(1'000'000);

        call_with_templates([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto rng = ankerl::nanobench::Rng{};

            auto vec = Vector{text};

            bench_rank.run(vector_name, [&]() {
                auto v = vec.rank(rng.bounded(text.size()));
                ankerl::nanobench::doNotOptimizeAway(v);
            });
        }, AllBitvectors{});
    }
}

TEST_CASE("benchmark bit vectors memory consumption", "[bitvector][size]") {
    BenchSize benchSize;
    benchSize.baseSize = 1.;

    SECTION("benchmarking") {
        call_with_templates([&]<typename Vector>() {

            auto vector_name = getName<Vector>();
            INFO(vector_name);

            auto& text = generateText();

            auto vec = Vector{text};
            if constexpr (requires() { vec.space_usage(); }) {
                auto s = vec.space_usage();
                benchSize.addEntry({
                    .name = vector_name,
                    .size = s,
                    .text_size = text.size(),
                    .bits_per_char = (s*8)/double(text.size())
                });
            } else {
                auto ofs     = std::stringstream{};
                auto archive = cereal::BinaryOutputArchive{ofs};
                archive(vec);
                auto s = ofs.str().size();
                benchSize.addEntry({
                    .name = vector_name,
                    .size = s,
                    .text_size = text.size(),
                    .bits_per_char = (s*8)/double(text.size())
                });
            }
        }, AllBitvectors{});
    }
}
