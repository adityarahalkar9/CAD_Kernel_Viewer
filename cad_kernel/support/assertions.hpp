// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : support/Assertions.hpp
//  Purpose : CAD_ASSERT / CAD_VERIFY / CAD_UNREACHABLE
//
//  CAD_ASSERT      internal invariants - compiled out with NDEBUG
//  CAD_VERIFY      always-on checks for untrusted input (loaded files, ...)
//  CAD_UNREACHABLE marks impossible control flow - aborts if ever reached
//
//  Both macros accept an optional string message:
//      CAD_ASSERT(index < size, "index out of range");
// ============================================================================
#pragma once

#include <cstdio>
#include <cstdlib>
#include <string_view>

// ---------------------------------------------------------------------------
// Debug break
// ---------------------------------------------------------------------------
#if defined(_MSC_VER)
#define CAD_DEBUGBREAK() __debugbreak()
#elif defined(__has_builtin)
#if __has_builtin(__builtin_trap)
#define CAD_DEBUGBREAK() __builtin_trap()
#else
#define CAD_DEBUGBREAK() ::std::abort()
#endif
#else
#define CAD_DEBUGBREAK() ::std::abort()
#endif

// ---------------------------------------------------------------------------
// Enable switch
// ---------------------------------------------------------------------------
#if !defined(NDEBUG) && !defined(CAD_DISABLE_ASSERTS)
#define CAD_ASSERT_ENABLED 1
#else
#define CAD_ASSERT_ENABLED 0
#endif

namespace cad::detail {

    // Prints a diagnostic, then terminates. In debug builds the debug break fires
    // first so an attached debugger stops exactly at the failure site; continuing
    // execution lands on abort().
    [[noreturn]] inline void report_and_abort(
        std::string_view kind,
        std::string_view expression,
        std::string_view message,
        const char* file,
        int               line)
    {
        std::fprintf(stderr,
                     "\n%s failed at %s:%d\n"
                     "  expression: %s\n",
                     kind.data(), file, line, expression.data());

        if(!message.empty()) {
            std::fprintf(stderr, "  message   : %s\n", message.data());
        }

        std::fflush(stderr);

#if CAD_ASSERT_ENABLED
        CAD_DEBUGBREAK();
#endif

        std::abort();
    }

} // namespace cad::detail

// ---------------------------------------------------------------------------
// Public macros
// ---------------------------------------------------------------------------
#if CAD_ASSERT_ENABLED
#define CAD_ASSERT(expression, ...)                                        \
        do {                                                                   \
            if (static_cast<bool>(expression) == false) {                      \
                ::cad::detail::report_and_abort(                               \
                    "CAD_ASSERT", #expression,                                 \
                    ::std::string_view("" __VA_OPT__(,) __VA_ARGS__),          \
                    __FILE__, __LINE__);                                       \
            }                                                                  \
        } while (false)
#else
#define CAD_ASSERT(expression, ...) static_cast<void>(0)
#endif

// Always-on validation for data that can come from outside the program
// (mesh files, user input). Never compiled out.
#define CAD_VERIFY(expression, ...)                                            \
    do {                                                                       \
        if (static_cast<bool>(expression) == false) {                          \
            ::cad::detail::report_and_abort(                                   \
                "CAD_VERIFY", #expression,                                     \
                ::std::string_view("" __VA_OPT__(,) __VA_ARGS__),              \
                __FILE__, __LINE__);                                           \
        }                                                                      \
    } while (false)

#define CAD_UNREACHABLE(...)                                                   \
    ::cad::detail::report_and_abort(                                           \
        "CAD_UNREACHABLE", "unreachable code path reached",                    \
        ::std::string_view("" __VA_OPT__(,) __VA_ARGS__),                      \
        __FILE__, __LINE__)