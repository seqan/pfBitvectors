# SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(pfBitvectors_test_utils)

add_library(${PROJECT_NAME} INTERFACE)
add_library(pfBitvectors::test_utils ALIAS ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} INTERFACE
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
