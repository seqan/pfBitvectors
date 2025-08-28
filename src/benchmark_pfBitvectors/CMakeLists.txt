# SPDX-FileCopyrightText: 2023 Gottlieb+Freitag <info@gottliebtfreitag.de>
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(benchmark_pfBitvectors)

add_executable(${PROJECT_NAME}
    benchmark_bitvectors.cpp
    benchmark_strings_alphabet_4.cpp
    benchmark_strings_alphabet_5.cpp
    benchmark_strings_alphabet_16.cpp
    benchmark_strings_alphabet_21.cpp
    benchmark_strings_alphabet_255.cpp
#    benchmark_strings_alphabet_4096.cpp
#    benchmark_strings_alphabet_16384.cpp
#    benchmark_strings_alphabet_65536.cpp
)
target_link_libraries(${PROJECT_NAME} PUBLIC
    Catch2::Catch2WithMain
    fmt::fmt
    cereal::cereal
    pfBitvectors::pfBitvectors
    pfBitvectors::externalLibsAdapter
    pfBitvectors::test_utils
    nanobench::nanobench
    reflect::reflect
)

if (WIN32)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/bigobj /EHsc")
elseif (APPLE)
else ()
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-Wall -Werror -Wpedantic")
endif ()


set_target_properties(${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
