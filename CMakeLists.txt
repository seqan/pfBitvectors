# SPDX-FileCopyrightText: 2023 Gottlieb+Freitag <info@gottliebtfreitag.de>
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(pfBitvectors LANGUAGES CXX C
            DESCRIPTION "pfBitvectors -- bit vectors and strings with rank support")

option(PFBITVECTORS_BENCHMARKS     "build benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(PFBITVECTORS_TESTS          "build tests"      ${PROJECT_IS_TOP_LEVEL})
option(PFBITVECTORS_USE_SDSL       "Activate occ tables that use SDSL as a backend" ${PROJECT_IS_TOP_LEVEL})
option(PFBITVECTORS_USE_PASTA      "Include PASTA (required for bitvector benchmarks)" ${PROJECT_IS_TOP_LEVEL})
option(PFBITVECTORS_USE_SUX        "Include SUX (required for bitvector benchmarks)" ${PROJECT_IS_TOP_LEVEL})
option(PFBITVECTORS_USE_AWFMINDEX  "Include AWFMIndex (only used in string benchmark)" ${PROJECT_IS_TOP_LEVEL})

if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(PFBITVECTORS_USE_SDSL OFF)
    set(PFBITVECTORS_USE_PASTA OFF)
    set(PFBITVECTORS_USE_SUX OFF)
    set(PFBITVECTORS_USE_AWFMINDEX OFF)
elseif (APPLE)
    set(PFBITVECTORS_USE_SDSL OFF)
    set(PFBITVECTORS_USE_PASTA OFF)
    set(PFBITVECTORS_USE_SUX OFF)
    set(PFBITVECTORS_USE_AWFMINDEX OFF)
elseif (EMSCRIPTEN)
    set(PFBITVECTORS_USE_SDSL OFF)
    set(PFBITVECTORS_USE_PASTA OFF)
    set(PFBITVECTORS_USE_SUX OFF)
    set(PFBITVECTORS_USE_AWFMINDEX OFF)
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(PFBITVECTORS_USE_SDSL OFF)
endif ()

include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/cpm.dependencies@1.0.1")
CPMLoadDependenciesFile("${CMAKE_CURRENT_SOURCE_DIR}/cpm.dependencies")

add_subdirectory(src/pfBitvectors)
add_subdirectory(src/pfBitvectors_externalLibsAdapter)

if (PFBITVECTORS_BENCHMARKS OR PFBITVECTORS_TESTS)
    add_subdirectory(src/pfBitvectors_test_utils)
endif ()

if (PFBITVECTORS_BENCHMARKS)
    add_subdirectory(src/benchmark_pfBitvectors)
endif ()
if (PFBITVECTORS_TESTS)
    enable_testing()
    add_subdirectory(src/test_pfBitvectors)
    add_subdirectory(src/test_header)
endif ()

