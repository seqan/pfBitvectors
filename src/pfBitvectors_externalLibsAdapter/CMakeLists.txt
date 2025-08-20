# SPDX-FileCopyrightText: 2023 Gottlieb+Freitag <info@gottliebtfreitag.de>
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(pfBitvectors_externalLibsAdapter)

add_library(${PROJECT_NAME} INTERFACE)
add_library(pfBitvectors::externalLibsAdapter ALIAS ${PROJECT_NAME})
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)

target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/..>
    $<INSTALL_INTERFACE:include>
)

if (PFBITVECTORS_USE_PASTA)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DPFBITVECTORS_USE_PASTA)
    target_link_libraries(${PROJECT_NAME} INTERFACE pasta_bit_vector)
endif ()

if (PFBITVECTORS_USE_SUX)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DPFBITVECTORS_USE_SUX)
    target_link_libraries(${PROJECT_NAME} INTERFACE sux::sux)
endif ()

if (PFBITVECTORS_USE_SDSL)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DPFBITVECTORS_USE_SDSL)
    target_link_libraries(${PROJECT_NAME} INTERFACE sdsl-lite::sdsl-lite)
endif ()

if (PFBITVECTORS_USE_AWFMINDEX)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DPFBITVECTORS_USE_AWFMINDEX)
    target_link_libraries(${PROJECT_NAME} INTERFACE AwFmIndex::AwFmIndex)
    target_compile_options(${PROJECT_NAME} INTERFACE -mavx2)
endif()

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
