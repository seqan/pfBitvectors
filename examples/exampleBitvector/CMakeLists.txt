
# SPDX-FileCopyrightText: 2023 Gottlieb+Freitag <info@gottliebtfreitag.de>
# SPDX-License-Identifier: CC0-1.0
cmake_minimum_required (VERSION 3.25)

project(example LANGUAGES CXX DESCRIPTION "example of SeqAns bitvectors")

include(cmake/CPM.cmake)
CPMAddPackage("gh:seqan/pfBitvectors#main")

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} seqan::pfb)
