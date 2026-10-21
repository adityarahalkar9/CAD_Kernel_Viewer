// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : math/Vec4.hpp
//  Purpose : 4-component vector (homogeneous coordinates, RGBA color)
//
//  CAD semantics helpers: Vec4::point(w = 1) and Vec4::direction(w = 0)
//  encode point/vector distinction used by every transform in Part 3.
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
#include "Vec2.hpp"
#include "Vec3.hpp"

namespace cad {

    template <typename T>
    struct Vec4 {
        static_assert(std::is_arithmetic_v<T>, "Vec4<T> requires an arithmetic T");

        using value_type = T;
        static constexpr int kDimension = 4;

        // -- data -----------------------------------------------------------------
        T x{};
        T y{};
        T z{};
        T w{};

        // -- construction -----------------------------------------------------------
        constexpr Vec4() noexcept = default;

        constexpr Vec4(T ix, T iy, T iz, T iw) noexcept
            : x(ix), y(iy), z(iz), w(iw) {
        }

        template <typename U>
            requires (!std::is_same_v<U, T>)
        explicit(!std::is_convertible_v<U, T>)
            constexpr Vec4(const Vec4<U>& other) noexcept
            : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)),
            z(static_cast<T>(other.z)), w(static_cast<T>(other.w)) {
        }

        // Dimension lift + w is always explicit.
        constexpr explicit Vec4(const Vec3<T>& v, T iw) noexcept
            : x(v.x), y(v.y), z(v.z), w(iw) {
        }

        // Homogeneous semantics: points carry w = 1, directions carry w = 0.
        static constexpr Vec4 point(const Vec3<T>& p) noexcept {
            return Vec4(p, T(1));
        }

        static constexpr Vec4 direction(const Vec3<T>& d) noexcept {
            return Vec4(d, T(0));
        }

        static constexpr Vec4 Zero() noexcept  { return Vec4{}; }
        static constexpr Vec4 One() noexcept   { return Vec4(T(1), T(1), T(1), T(1)); }
        static constexpr Vec4 UnitX() noexcept { return Vec4(T(1), T(0), T(0), T(0)); }
        static constexpr Vec4 UnitY() noexcept { return Vec4(T(0), T(1), T(0), T(0)); }
        static constexpr Vec4 UnitZ() noexcept { return Vec4(T(0), T(0), T(1), T(0)); }
        static constexpr Vec4 UnitW() noexcept { return Vec4(T(0), T(0), T(0), T(1)); }

        // -- access ------------------------------------------------------------------
        [[nodiscard]] constexpr T& operator[](int index) noexcept {
            CAD_ASSERT(index >= 0 && index < 4, "Vec4 index out of range");
            return index == 0 ? x : (index == 1 ? y : (index == 2 ? z : w));
        }

        [[nodiscard]] constexpr const T& operator[](int index) const noexcept {
            CAD_ASSERT(index >= 0 && index < 4, "Vec4 index out of range");
            return index == 0 ? x : (index == 1 ? y : (index == 2 ? z : w));
        }

        [[nodiscard]] constexpr Vec2<T> xy() const noexcept { return Vec2<T>(x, y); }
        [[nodiscard]] constexpr Vec3<T> xyz() const noexcept { return Vec3<T>(x, y, z); }

        // -- compound arithmetic --------------------------------------------------------
        constexpr Vec4& operator+=(const Vec4& r) noexcept { x += r.x; y += r.y; z += r.z; w += r.w; return *this; }
        constexpr Vec4& operator-=(const Vec4& r) noexcept { x -= r.x; y -= r.y; z -= r.z; w -= r.w; return *this; }
        constexpr Vec4& operator*=(T s) noexcept { x *= s; y *= s; z *= s; w *= s; return *this; }

        constexpr Vec4& operator/=(T s) noexcept {
            CAD_ASSERT(s != T(0), "Vec4 division by zero");
            x /= s;
            y /= s;
            z /= s;
            w /= s;
            return *this;
        }

        // -- unary -----------------------------------------------------------------
        [[nodiscard]] constexpr Vec4 operator-() const noexcept { return Vec4(-x, -y, -z, -w); }

        // -- products & norms ---------------------------------------------------------
        [[nodiscard]] constexpr T dot(const Vec4& r) const noexcept {
            return x * r.x + y * r.y + z * r.z + w * r.w;
        }

        [[nodiscard]] constexpr T lengthSquared() const noexcept {
            return x * x + y * y + z * z + w * w;
        }

        [[nodiscard]] constexpr T length() const noexcept
            requires std::is_floating_point_v<T>
        {
            return cad::sqrt(lengthSquared());
        }

        [[nodiscard]] constexpr bool isZero() const noexcept {
            return x == T(0) && y == T(0) && z == T(0) && w == T(0);
        }

        [[nodiscard]] constexpr Vec4 normalized() const noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(lengthSquared() > T(0), "Vec4::normalized on zero-length vector");
            const T inv = T(1) / cad::sqrt(lengthSquared());
            return Vec4(x * inv, y * inv, z * inv, w * inv);
        }

        [[nodiscard]] constexpr Vec4 normalizedOr(const Vec4& fallback) const noexcept
            requires std::is_floating_point_v<T>
        {
            const T len2 = lengthSquared();
            if(len2 == T(0)) {
                return fallback;
            }
            const T inv = T(1) / cad::sqrt(len2);
            return Vec4(x * inv, y * inv, z * inv, w * inv);
        }

        [[nodiscard]] constexpr bool isNormalized(
            T tol = GeometryTolerance<T>::unit_deviation) const noexcept
            requires std::is_floating_point_v<T>
        {
            return almostZero(lengthSquared() - T(1), tol);
        }

        // Homogeneous divide: (x/w, y/w, z/w).
        [[nodiscard]] constexpr Vec3<T> perspectiveDivide() const noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(w != T(0), "perspectiveDivide: w == 0");
            const T inv = T(1) / w;
            return Vec3<T>(x * inv, y * inv, z * inv);
        }

        // -- component reducers -----------------------------------------------------
        [[nodiscard]] constexpr T minComponent() const noexcept {
            const T m12 = y < z ? y : z;
            const T m = x < m12 ? x : m12;
            return m < w ? m : w;
        }

        [[nodiscard]] constexpr T maxComponent() const noexcept {
            const T m12 = y > z ? y : z;
            const T m = x > m12 ? x : m12;
            return m > w ? m : w;
        }

        [[nodiscard]] constexpr int minAxis() const noexcept {
            const T m = minComponent();
            if(m == x) return 0;
            if(m == y) return 1;
            if(m == z) return 2;
            return 3;
        }

        [[nodiscard]] constexpr int maxAxis() const noexcept {
            const T m = maxComponent();
            if(m == x) return 0;
            if(m == y) return 1;
            if(m == z) return 2;
            return 3;
        }

        // -- componentwise maps --------------------------------------------------------
        [[nodiscard]] constexpr Vec4 cwAbs() const noexcept {
            return Vec4(x < T(0) ? -x : x,
                        y < T(0) ? -y : y,
                        z < T(0) ? -z : z,
                        w < T(0) ? -w : w);
        }

        [[nodiscard]] constexpr Vec4 cwMin(const Vec4& r) const noexcept {
            return Vec4(x < r.x ? x : r.x,
                        y < r.y ? y : r.y,
                        z < r.z ? z : r.z,
                        w < r.w ? w : r.w);
        }

        [[nodiscard]] constexpr Vec4 cwMax(const Vec4& r) const noexcept {
            return Vec4(x > r.x ? x : r.x,
                        y > r.y ? y : r.y,
                        z > r.z ? z : r.z,
                        w > r.w ? w : r.w);
        }

        [[nodiscard]] Vec4 cwFloor() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec4(std::floor(x), std::floor(y), std::floor(z), std::floor(w));
        }

        [[nodiscard]] Vec4 cwCeil() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec4(std::ceil(x), std::ceil(y), std::ceil(z), std::ceil(w));
        }

        [[nodiscard]] Vec4 cwRound() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec4(std::round(x), std::round(y), std::round(z), std::round(w));
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return std::isfinite(x) && std::isfinite(y)
                && std::isfinite(z) && std::isfinite(w);
        }
    };

    // -- type aliases ------------------------------------------------------------
    using Vec4f = Vec4<float>;
    using Vec4d = Vec4<double>;
    using Vec4i = Vec4<std::int32_t>;

    // -- free operators (scalar ops promote to the wider type) ----------------------
    template <typename T>
    [[nodiscard]] constexpr Vec4<T> operator+(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        return Vec4<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec4<T> operator-(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        return Vec4<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec4<std::common_type_t<T, U>> operator*(const Vec4<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        return Vec4<R>(R(v.x) * R(s), R(v.y) * R(s), R(v.z) * R(s), R(v.w) * R(s));
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec4<std::common_type_t<T, U>> operator*(U s, const Vec4<T>& v) noexcept {
        return v * s;
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec4<std::common_type_t<T, U>> operator/(const Vec4<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        CAD_ASSERT(R(s) != R(0), "Vec4 division by zero");
        return Vec4<R>(R(v.x) / R(s), R(v.y) / R(s), R(v.z) / R(s), R(v.w) / R(s));
    }

    template <typename T>
    [[nodiscard]] constexpr Vec4<T> cwMul(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        return Vec4<T>(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec4<T> cwDiv(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        CAD_ASSERT(b.x != T(0) && b.y != T(0) && b.z != T(0) && b.w != T(0),
                   "cwDiv: zero component");
        return Vec4<T>(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
    }

    // -- exact comparison + lexicographic ordering -----------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator<(const Vec4<T>& a, const Vec4<T>& b) noexcept {
        if(a.x != b.x) return a.x < b.x;
        if(a.y != b.y) return a.y < b.y;
        if(a.z != b.z) return a.z < b.z;
        return a.w < b.w;
    }

    // -- free-function spellings of the hot members -------------------------------------
    template <typename T>
    [[nodiscard]] constexpr T dot(const Vec4<T>& a, const Vec4<T>& b) noexcept { return a.dot(b); }

    template <typename T>
    [[nodiscard]] constexpr T lengthSquared(const Vec4<T>& v) noexcept { return v.lengthSquared(); }

    template <typename T>
    [[nodiscard]] constexpr T length(const Vec4<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.length();
    }

    template <typename T>
    [[nodiscard]] constexpr Vec4<T> normalize(const Vec4<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.normalized();
    }

    template <typename T>
    [[nodiscard]] constexpr Vec4<T> lerp(const Vec4<T>& a, const Vec4<T>& b, T t) noexcept
        requires std::is_floating_point_v<T>
    {
        return Vec4<T>(cad::lerp(a.x, b.x, t),
                       cad::lerp(a.y, b.y, t),
                       cad::lerp(a.z, b.z, t),
                       cad::lerp(a.w, b.w, t));
    }

    // -- approximate equality ---------------------------------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Vec4<T>& a, const Vec4<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.x, b.x, relTol, absTol)
            && almostEqual(a.y, b.y, relTol, absTol)
            && almostEqual(a.z, b.z, relTol, absTol)
            && almostEqual(a.w, b.w, relTol, absTol);
    }

    // -- stream output -------------------------------------------------------------------
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Vec4<T>& v) {
        return os << std::format("({}, {}, {}, {})", v.x, v.y, v.z, v.w);
    }

} // namespace cad

// -- std::formatter / std::hash specializations --------------------------------
namespace std {

    template <typename T>
    struct formatter<cad::Vec4<T>, char> {
        cad::detail::FormatSpec spec;

        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }

        auto format(const cad::Vec4<T>& value, std::format_context& ctx) const {
            return cad::detail::emitBracketed(ctx, "(", ")", spec.view(),
                                              value.x, value.y, value.z, value.w);
        }
    };

    template <typename T>
    struct hash<cad::Vec4<T>> {
        [[nodiscard]] std::size_t operator()(const cad::Vec4<T>& v) const noexcept {
            std::size_t seed = std::hash<T>{}(v.x);
            seed = cad::hashCombine(seed, std::hash<T>{}(v.y));
            seed = cad::hashCombine(seed, std::hash<T>{}(v.z));
            return cad::hashCombine(seed, std::hash<T>{}(v.w));
        }
    };

} // namespace std