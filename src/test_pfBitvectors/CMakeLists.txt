# SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(test_pfBitvectors)

add_executable(${PROJECT_NAME}
    test_bitvectors.cpp
    test_strings.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2WithMain
    pfBitvectors::pfBitvectors
    pfBitvectors::externalLibsAdapter
    pfBitvectors::test_utils
    cereal::cereal
    reflect::reflect
)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

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
