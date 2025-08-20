// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>
#include <variant>

#if __has_include(<cxxabi.h>)
#include <cxxabi.h>
namespace {
template <typename T>
auto getName() {
    int     status;
    auto realname = abi::__cxa_demangle(typeid(T).name(), nullptr, nullptr, &status);
    auto str = std::string{realname};
    std::free(realname);
    return str;
}
}
#else
#include <reflect>
namespace {
template <typename T>
auto getName() {
    return std::string{reflect::type_name<T>()};
}
}
#endif

namespace {

template <typename... T>
void call_with_templates(auto f, std::variant<T...> const&) {
    auto g = [&]<typename T2>() {
        if constexpr (!std::same_as<T2, std::monostate>) {
            f.template operator()<T2>();
        }
    };
    (g.template operator()<T>(), ...);
}

template <template <size_t> class... T>
struct Variant {};

template <size_t>
struct Delimiter {};

template <template <size_t> class... T>
void call_with_templates(auto f, Variant<T...> const&) {
    auto g = [&]<template <size_t> typename T2>() {
        if constexpr (!std::same_as<T2<4>, Delimiter<4>>) {
            f.template operator()<T2>();
        }
    };
    (g.template operator()<T>(), ...);
}

template <typename ...Ts>
struct AppendImpl;

template <typename... T1, typename... Ts>
struct AppendImpl<std::variant<T1...>, Ts...> {
    using type = std::variant<T1..., Ts...>;
};


template <typename T1, typename ...Ts>
using Append = AppendImpl<T1, Ts...>::type;

}
