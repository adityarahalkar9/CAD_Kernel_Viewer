// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : math/Vec2.hpp
//  Purpose : 2-component vector - constexpr-first, std::format/std::hash aware
//
//  Also hosts the cad::detail formatting plumbing shared by Vec2/Vec3/Vec4
//  (this is the root of the vector include chain).
//
//  Design notes:
//   * default constructor zero-initializes (deterministic, CAD-friendly)
//   * exact == / != (approximate comparison lives in almostEqual)
//   * scalar * and / PROMOTE to the wider type - no silent precision loss
//   * length()/normalized() are constexpr via cad::sqrt
//   * floating-only members are constrained out for integer T
// ============================================================================
#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "Constants.hpp"
#include "Utils.hpp"
#include "../support/Assertions.hpp"

namespace cad {

    template <typename T>
    struct Vec2 {
        static_assert(std::is_arithmetic_v<T>, "Vec2<T> requires an arithmetic T");

        using value_type = T;
        static constexpr int kDimension = 2;

        // -- data -----------------------------------------------------------------
        T x{};
        T y{};

        // -- construction -----------------------------------------------------------
        constexpr Vec2() noexcept = default;

        constexpr Vec2(T ix, T iy) noexcept
            : x(ix), y(iy) {
        }

        // Cross-type conversion: implicit when lossless (float -> double),
        // explicit when narrowing (double -> float).
        template <typename U>
            requires (!std::is_same_v<U, T>)
        explicit(!std::is_convertible_v<U, T>)
            constexpr Vec2(const Vec2<U>& other) noexcept
            : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)) {
        }

        static constexpr Vec2 Zero() noexcept  { return Vec2{}; }
        static constexpr Vec2 One() noexcept   { return Vec2(T(1), T(1)); }
        static constexpr Vec2 UnitX() noexcept { return Vec2(T(1), T(0)); }
        static constexpr Vec2 UnitY() noexcept { return Vec2(T(0), T(1)); }

        // -- access ------------------------------------------------------------------
        [[nodiscard]] constexpr T& operator[](int index) noexcept {
            CAD_ASSERT(index == 0 || index == 1, "Vec2 index out of range");
            return index == 0 ? x : y;
        }

        [[nodiscard]] constexpr const T& operator[](int index) const noexcept {
            CAD_ASSERT(index == 0 || index == 1, "Vec2 index out of range");
            return index == 0 ? x : y;
        }

        // -- compound arithmetic --------------------------------------------------------
        constexpr Vec2& operator+=(const Vec2& r) noexcept { x += r.x; y += r.y; return *this; }
        constexpr Vec2& operator-=(const Vec2& r) noexcept { x -= r.x; y -= r.y; return *this; }
        constexpr Vec2& operator*=(T s) noexcept { x *= s; y *= s; return *this; }

        constexpr Vec2& operator/=(T s) noexcept {
            CAD_ASSERT(s != T(0), "Vec2 division by zero");
            x /= s;
            y /= s;
            return *this;
        }

        // -- unary -----------------------------------------------------------------
        [[nodiscard]] constexpr Vec2 operator-() const noexcept { return Vec2(-x, -y); }

        // -- products & norms ---------------------------------------------------------
        [[nodiscard]] constexpr T dot(const Vec2& r) const noexcept { return x * r.x + y * r.y; }

        // 2D pseudo cross product (the z of the 3D cross of the embedded vectors):
        //   > 0 -> r is CCW of *this, < 0 -> CW, == 0 -> parallel.
        [[nodiscard]] constexpr T cross(const Vec2& r) const noexcept { return x * r.y - y * r.x; }

        // Counterclockwise quarter turn.
        [[nodiscard]] constexpr Vec2 perp() const noexcept { return Vec2(-y, x); }

        [[nodiscard]] constexpr T lengthSquared() const noexcept { return x * x + y * y; }

        [[nodiscard]] constexpr T length() const noexcept
            requires std::is_floating_point_v<T>
        {
            return cad::sqrt(lengthSquared());
        }

        [[nodiscard]] constexpr bool isZero() const noexcept {
            return x == T(0) && y == T(0);
        }

        [[nodiscard]] constexpr Vec2 normalized() const noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(lengthSquared() > T(0), "Vec2::normalized on zero-length vector");
            const T inv = T(1) / cad::sqrt(lengthSquared());
            return Vec2(x * inv, y * inv);
        }

        [[nodiscard]] constexpr Vec2 normalizedOr(const Vec2& fallback) const noexcept
            requires std::is_floating_point_v<T>
        {
            const T len2 = lengthSquared();
            if(len2 == T(0)) {
                return fallback;
            }
            const T inv = T(1) / cad::sqrt(len2);
            return Vec2(x * inv, y * inv);
        }

        [[nodiscard]] constexpr bool isNormalized(
            T tol = GeometryTolerance<T>::unit_deviation) const noexcept
            requires std::is_floating_point_v<T>
        {
            return almostZero(lengthSquared() - T(1), tol);
        }

        // -- component reducers -----------------------------------------------------
        [[nodiscard]] constexpr T minComponent() const noexcept { return x < y ? x : y; }
        [[nodiscard]] constexpr T maxComponent() const noexcept { return x > y ? x : y; }
        [[nodiscard]] constexpr int minAxis() const noexcept { return x < y ? 0 : 1; }
        [[nodiscard]] constexpr int maxAxis() const noexcept { return x > y ? 0 : 1; }

        // -- componentwise maps (all constexpr) ----------------------------------------
        [[nodiscard]] constexpr Vec2 cwAbs() const noexcept {
            return Vec2(x < T(0) ? -x : x, y < T(0) ? -y : y);
        }
        [[nodiscard]] constexpr Vec2 cwMin(const Vec2& r) const noexcept {
            return Vec2(x < r.x ? x : r.x, y < r.y ? y : r.y);
        }
        [[nodiscard]] constexpr Vec2 cwMax(const Vec2& r) const noexcept {
            return Vec2(x > r.x ? x : r.x, y > r.y ? y : r.y);
        }

        // -- floating-only runtime maps --------------------------------------------------
        [[nodiscard]] Vec2 cwFloor() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec2(std::floor(x), std::floor(y));
        }

        [[nodiscard]] Vec2 cwCeil() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec2(std::ceil(x), std::ceil(y));
        }

        [[nodiscard]] Vec2 cwRound() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec2(std::round(x), std::round(y));
        }

        // -- directions & angles (floating only) -------------------------------------------
        // Signed angle from *this to r, in (-pi, pi]. Robust for non-unit inputs.
        [[nodiscard]] T angleTo(const Vec2& r) const noexcept
            requires std::is_floating_point_v<T>
        {
            return std::atan2(cross(r), dot(r));
        }

        [[nodiscard]] Vec2 rotated(T angleRadians) const noexcept
            requires std::is_floating_point_v<T>
        {
            const T c = std::cos(angleRadians);
            const T s = std::sin(angleRadians);
            return Vec2(x * c - y * s, x * s + y * c);
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return std::isfinite(x) && std::isfinite(y);
        }
    };

    // -- type aliases ------------------------------------------------------------
    using Vec2f = Vec2<float>;
    using Vec2d = Vec2<double>;
    using Vec2i = Vec2<std::int32_t>;

    // -- free operators (scalar ops promote to the wider type) ----------------------
    template <typename T>
    [[nodiscard]] constexpr Vec2<T> operator+(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        return Vec2<T>(a.x + b.x, a.y + b.y);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec2<T> operator-(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        return Vec2<T>(a.x - b.x, a.y - b.y);
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec2<std::common_type_t<T, U>> operator*(const Vec2<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        return Vec2<R>(R(v.x) * R(s), R(v.y) * R(s));
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec2<std::common_type_t<T, U>> operator*(U s, const Vec2<T>& v) noexcept {
        return v * s;
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec2<std::common_type_t<T, U>> operator/(const Vec2<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        CAD_ASSERT(R(s) != R(0), "Vec2 division by zero");
        return Vec2<R>(R(v.x) / R(s), R(v.y) / R(s));
    }

    // Componentwise multiply/divide (no operator* between vectors - on purpose:
    // dot/cross are explicit, so v * w can never silently mean the wrong thing).
    template <typename T>
    [[nodiscard]] constexpr Vec2<T> cwMul(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        return Vec2<T>(a.x * b.x, a.y * b.y);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec2<T> cwDiv(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        CAD_ASSERT(b.x != T(0) && b.y != T(0), "cwDiv: zero component");
        return Vec2<T>(a.x / b.x, a.y / b.y);
    }

    // -- exact comparison + lexicographic ordering -----------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        return a.x == b.x && a.y == b.y;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator<(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        if(a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    }

    // -- free-function spellings of the hot members -------------------------------------
    template <typename T>
    [[nodiscard]] constexpr T dot(const Vec2<T>& a, const Vec2<T>& b) noexcept { return a.dot(b); }

    template <typename T>
    [[nodiscard]] constexpr T cross(const Vec2<T>& a, const Vec2<T>& b) noexcept { return a.cross(b); }

    template <typename T>
    [[nodiscard]] constexpr T lengthSquared(const Vec2<T>& v) noexcept { return v.lengthSquared(); }

    template <typename T>
    [[nodiscard]] constexpr T length(const Vec2<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.length();
    }

    template <typename T>
    [[nodiscard]] constexpr Vec2<T> normalize(const Vec2<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.normalized();
    }

    template <typename T>
    [[nodiscard]] constexpr T distanceSquared(const Vec2<T>& a, const Vec2<T>& b) noexcept {
        const Vec2<T> d = b - a;
        return d.lengthSquared();
    }

    template <typename T>
    [[nodiscard]] constexpr T distance(const Vec2<T>& a, const Vec2<T>& b) noexcept
        requires std::is_floating_point_v<T>
    {
        return cad::sqrt(distanceSquared(a, b));
    }

    template <typename T>
    [[nodiscard]] constexpr Vec2<T> lerp(const Vec2<T>& a, const Vec2<T>& b, T t) noexcept
        requires std::is_floating_point_v<T>
    {
        return Vec2<T>(cad::lerp(a.x, b.x, t), cad::lerp(a.y, b.y, t));
    }

    // -- approximate equality (componentwise) -----------------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Vec2<T>& a, const Vec2<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.x, b.x, relTol, absTol)
            && almostEqual(a.y, b.y, relTol, absTol);
    }

    // -- stream output (also what GoogleTest prints on failures) ------------------------
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
        return os << std::format("({}, {})", v.x, v.y);
    }

    // ---------------------------------------------------------------------------
    // std::format plumbing shared by Vec2/Vec3/Vec4
    // ---------------------------------------------------------------------------
    namespace detail {

        // Compile-time buffer holding a std::format spec (the part between ':' and
        // '}'). POD on purpose: formatter::parse stays constexpr-friendly with zero
        // allocations on every toolchain.
        struct FormatSpec {
            static constexpr std::size_t kMax = 24;

            char         chars[kMax]{};
            std::size_t  count = 0;

            constexpr void append(char c) {
                CAD_ASSERT(count < kMax, "format spec too long");
                chars[count++] = c;
            }

            [[nodiscard]] std::string_view view() const noexcept {
                return std::string_view(chars, count);
            }
        };

        // Consumes spec characters up to (not including) '}' and returns an iterator
        // positioned at the '}'. Throws std::format_error on a nested '{' or a
        // missing '}' - which surfaces as a compile error inside std::format("...").
        template <typename ParseContext>
        constexpr auto parseSpec(FormatSpec& spec, ParseContext& ctx) {
            auto it = ctx.begin();
            auto end = ctx.end();

            while(it != end && *it != '}') {
                if(*it == '{') {
                    throw std::format_error("nested '{' in Vec format spec");
                }
                spec.append(*it);
                ++it;
            }
            if(it == end) {
                throw std::format_error("missing '}' in Vec format spec");
            }
            return it;
        }

        // Emits "open c1, c2, ... close" with the captured spec applied to every
        // component, e.g. "{:8.3f}" -> "(   1.000,    2.000)".
        template <typename... Components>
        auto emitBracketed(std::format_context& ctx,
                           std::string_view     open,
                           std::string_view     close,
                           std::string_view     spec,
                           const Components&... components)
        {
            std::string piece;
            piece.reserve(spec.size() + 3);
            piece += "{:";
            piece += spec;
            piece += '}';

            auto out = ctx.out();
            out = std::format_to(out, "{}", open);

            std::size_t index = 0;
            const auto  emitOne = [&](const auto& component) {
                if(index++ != 0) {
                    out = std::format_to(out, ", ");
                }
                out = std::vformat_to(out, piece, std::make_format_args(component));
                };
            (emitOne(components), ...);

            return std::format_to(out, "{}", close);
        }

    } // namespace detail

} // namespace cad

// -- std::formatter / std::hash specializations --------------------------------
namespace std {

    template <typename T>
    struct formatter<cad::Vec2<T>, char> {
        cad::detail::FormatSpec spec;

        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }

        auto format(const cad::Vec2<T>& value, std::format_context& ctx) const {
            return cad::detail::emitBracketed(ctx, "(", ")", spec.view(),
                                              value.x, value.y);
        }
    };

    template <typename T>
    struct hash<cad::Vec2<T>> {
        [[nodiscard]] std::size_t operator()(const cad::Vec2<T>& v) const noexcept {
            std::size_t seed = std::hash<T>{}(v.x);
            return cad::hashCombine(seed, std::hash<T>{}(v.y));
        }
    };

} // namespace std