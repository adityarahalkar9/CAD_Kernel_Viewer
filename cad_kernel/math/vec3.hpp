// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : math/Vec3.hpp
//  Purpose : 3-component vector - THE workhorse type of the kernel
//
//  Everything Vec2 offers plus: cross product, orthonormal basis construction
//  (Duff et al. 2017 - branchless and division-safe), angleBetween, reflect,
//  axis reducers used by the watertight ray-triangle kernel in Part 4.
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

namespace cad {

    template <typename T>
    struct Vec3 {
        static_assert(std::is_arithmetic_v<T>, "Vec3<T> requires an arithmetic T");

        using value_type = T;
        static constexpr int kDimension = 3;

        // -- data -----------------------------------------------------------------
        T x{};
        T y{};
        T z{};

        // -- construction -----------------------------------------------------------
        constexpr Vec3() noexcept = default;

        constexpr Vec3(T ix, T iy, T iz) noexcept
            : x(ix), y(iy), z(iz) {
        }

        template <typename U>
            requires (!std::is_same_v<U, T>)
        explicit(!std::is_convertible_v<U, T>)
            constexpr Vec3(const Vec3<U>& other) noexcept
            : x(static_cast<T>(other.x)), y(static_cast<T>(other.y)),
            z(static_cast<T>(other.z)) {
        }

        // Dimension lift from 2D is always explicit (z = 0 is a semantic choice).
        template <typename U>
        explicit constexpr Vec3(const Vec2<U>& v) noexcept
            : x(static_cast<T>(v.x)), y(static_cast<T>(v.y)), z(T(0)) {
        }

        static constexpr Vec3 Zero() noexcept  { return Vec3{}; }
        static constexpr Vec3 One() noexcept   { return Vec3(T(1), T(1), T(1)); }
        static constexpr Vec3 UnitX() noexcept { return Vec3(T(1), T(0), T(0)); }
        static constexpr Vec3 UnitY() noexcept { return Vec3(T(0), T(1), T(0)); }
        static constexpr Vec3 UnitZ() noexcept { return Vec3(T(0), T(0), T(1)); }

        // -- access ------------------------------------------------------------------
        [[nodiscard]] constexpr T& operator[](int index) noexcept {
            CAD_ASSERT(index >= 0 && index < 3, "Vec3 index out of range");
            return index == 0 ? x : (index == 1 ? y : z);
        }

        [[nodiscard]] constexpr const T& operator[](int index) const noexcept {
            CAD_ASSERT(index >= 0 && index < 3, "Vec3 index out of range");
            return index == 0 ? x : (index == 1 ? y : z);
        }

        [[nodiscard]] constexpr Vec2<T> xy() const noexcept { return Vec2<T>(x, y); }

        // -- compound arithmetic --------------------------------------------------------
        constexpr Vec3& operator+=(const Vec3& r) noexcept { x += r.x; y += r.y; z += r.z; return *this; }
        constexpr Vec3& operator-=(const Vec3& r) noexcept { x -= r.x; y -= r.y; z -= r.z; return *this; }
        constexpr Vec3& operator*=(T s) noexcept { x *= s; y *= s; z *= s; return *this; }

        constexpr Vec3& operator/=(T s) noexcept {
            CAD_ASSERT(s != T(0), "Vec3 division by zero");
            x /= s;
            y /= s;
            z /= s;
            return *this;
        }

        // -- unary -----------------------------------------------------------------
        [[nodiscard]] constexpr Vec3 operator-() const noexcept { return Vec3(-x, -y, -z); }

        // -- products & norms ---------------------------------------------------------
        [[nodiscard]] constexpr T dot(const Vec3& r) const noexcept {
            return x * r.x + y * r.y + z * r.z;
        }

        // Right-handed cross product. THE topology/normal workhorse.
        [[nodiscard]] constexpr Vec3 cross(const Vec3& r) const noexcept {
            return Vec3(y * r.z - z * r.y,
                        z * r.x - x * r.z,
                        x * r.y - y * r.x);
        }

        [[nodiscard]] constexpr T lengthSquared() const noexcept {
            return x * x + y * y + z * z;
        }

        [[nodiscard]] constexpr T length() const noexcept
            requires std::is_floating_point_v<T>
        {
            return cad::sqrt(lengthSquared());
        }

        [[nodiscard]] constexpr bool isZero() const noexcept {
            return x == T(0) && y == T(0) && z == T(0);
        }

        [[nodiscard]] constexpr Vec3 normalized() const noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(lengthSquared() > T(0), "Vec3::normalized on zero-length vector");
            const T inv = T(1) / cad::sqrt(lengthSquared());
            return Vec3(x * inv, y * inv, z * inv);
        }

        [[nodiscard]] constexpr Vec3 normalizedOr(const Vec3& fallback) const noexcept
            requires std::is_floating_point_v<T>
        {
            const T len2 = lengthSquared();
            if(len2 == T(0)) {
                return fallback;
            }
            const T inv = T(1) / cad::sqrt(len2);
            return Vec3(x * inv, y * inv, z * inv);
        }

        [[nodiscard]] constexpr bool isNormalized(
            T tol = GeometryTolerance<T>::unit_deviation) const noexcept
            requires std::is_floating_point_v<T>
        {
            return almostZero(lengthSquared() - T(1), tol);
        }

        // -- component reducers -----------------------------------------------------
        [[nodiscard]] constexpr T minComponent() const noexcept {
            return x < y ? (x < z ? x : z) : (y < z ? y : z);
        }

        [[nodiscard]] constexpr T maxComponent() const noexcept {
            return x > y ? (x > z ? x : z) : (y > z ? y : z);
        }

        [[nodiscard]] constexpr int minAxis() const noexcept {
            if(x <= y && x <= z) return 0;
            if(y <= z)           return 1;
            return 2;
        }

        [[nodiscard]] constexpr int maxAxis() const noexcept {
            if(x >= y && x >= z) return 0;
            if(y >= z)           return 1;
            return 2;
        }

        // -- componentwise maps --------------------------------------------------------
        [[nodiscard]] constexpr Vec3 cwAbs() const noexcept {
            return Vec3(x < T(0) ? -x : x,
                        y < T(0) ? -y : y,
                        z < T(0) ? -z : z);
        }

        [[nodiscard]] constexpr Vec3 cwMin(const Vec3& r) const noexcept {
            return Vec3(x < r.x ? x : r.x,
                        y < r.y ? y : r.y,
                        z < r.z ? z : r.z);
        }

        [[nodiscard]] constexpr Vec3 cwMax(const Vec3& r) const noexcept {
            return Vec3(x > r.x ? x : r.x,
                        y > r.y ? y : r.y,
                        z > r.z ? z : r.z);
        }

        [[nodiscard]] Vec3 cwFloor() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec3(std::floor(x), std::floor(y), std::floor(z));
        }

        [[nodiscard]] Vec3 cwCeil() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec3(std::ceil(x), std::ceil(y), std::ceil(z));
        }

        [[nodiscard]] Vec3 cwRound() const noexcept
            requires std::is_floating_point_v<T>
        {
            return Vec3(std::round(x), std::round(y), std::round(z));
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
        }
    };

    // -- type aliases ------------------------------------------------------------
    using Vec3f = Vec3<float>;
    using Vec3d = Vec3<double>;
    using Vec3i = Vec3<std::int32_t>;

    // -- free operators (scalar ops promote to the wider type) ----------------------
    template <typename T>
    [[nodiscard]] constexpr Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        return Vec3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        return Vec3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec3<std::common_type_t<T, U>> operator*(const Vec3<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        return Vec3<R>(R(v.x) * R(s), R(v.y) * R(s), R(v.z) * R(s));
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec3<std::common_type_t<T, U>> operator*(U s, const Vec3<T>& v) noexcept {
        return v * s;
    }

    template <typename T, typename U>
        requires std::is_arithmetic_v<U>
    [[nodiscard]] constexpr Vec3<std::common_type_t<T, U>> operator/(const Vec3<T>& v, U s) noexcept {
        using R = std::common_type_t<T, U>;
        CAD_ASSERT(R(s) != R(0), "Vec3 division by zero");
        return Vec3<R>(R(v.x) / R(s), R(v.y) / R(s), R(v.z) / R(s));
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> cwMul(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        return Vec3<T>(a.x * b.x, a.y * b.y, a.z * b.z);
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> cwDiv(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        CAD_ASSERT(b.x != T(0) && b.y != T(0) && b.z != T(0), "cwDiv: zero component");
        return Vec3<T>(a.x / b.x, a.y / b.y, a.z / b.z);
    }

    // -- exact comparison + lexicographic ordering -----------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator<(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        if(a.x != b.x) return a.x < b.x;
        if(a.y != b.y) return a.y < b.y;
        return a.z < b.z;
    }

    // -- free-function spellings of the hot members -------------------------------------
    template <typename T>
    [[nodiscard]] constexpr T dot(const Vec3<T>& a, const Vec3<T>& b) noexcept { return a.dot(b); }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b) noexcept { return a.cross(b); }

    template <typename T>
    [[nodiscard]] constexpr T lengthSquared(const Vec3<T>& v) noexcept { return v.lengthSquared(); }

    template <typename T>
    [[nodiscard]] constexpr T length(const Vec3<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.length();
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> normalize(const Vec3<T>& v) noexcept
        requires std::is_floating_point_v<T>
    {
        return v.normalized();
    }

    template <typename T>
    [[nodiscard]] constexpr T distanceSquared(const Vec3<T>& a, const Vec3<T>& b) noexcept {
        const Vec3<T> d = b - a;
        return d.lengthSquared();
    }

    template <typename T>
    [[nodiscard]] constexpr T distance(const Vec3<T>& a, const Vec3<T>& b) noexcept
        requires std::is_floating_point_v<T>
    {
        return cad::sqrt(distanceSquared(a, b));
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> lerp(const Vec3<T>& a, const Vec3<T>& b, T t) noexcept
        requires std::is_floating_point_v<T>
    {
        return Vec3<T>(cad::lerp(a.x, b.x, t),
                       cad::lerp(a.y, b.y, t),
                       cad::lerp(a.z, b.z, t));
    }

    // Unclamped angle in [0, pi] - robust for non-normalized inputs, unlike
    // acos(dot) which loses all precision near 0 and pi.
    template <typename T>
    [[nodiscard]] T angleBetween(const Vec3<T>& a, const Vec3<T>& b) noexcept
        requires std::is_floating_point_v<T>
    {
        return std::atan2(a.cross(b).length(), a.dot(b));
    }

    // Reflect incident v about a NORMALIZED normal n.
    template <typename T>
    [[nodiscard]] constexpr Vec3<T> reflect(const Vec3<T>& v, const Vec3<T>& n) noexcept
        requires std::is_floating_point_v<T>
    {
        return v - n * (T(2) * v.dot(n));
    }

    // -- approximate equality ---------------------------------------------------------
    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Vec3<T>& a, const Vec3<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.x, b.x, relTol, absTol)
            && almostEqual(a.y, b.y, relTol, absTol)
            && almostEqual(a.z, b.z, relTol, absTol);
    }

    // -- orthonormal basis -------------------------------------------------------------
    // "Building an Orthonormal Basis, Revisited" - Duff et al., JCGT 2017.
    // Branchless, no division-by-zero for any NORMALIZED n (the denominator
    // sign + n.z is never zero: it is 0 only when |n.z| > 1). Returns a
    // right-handed frame (tangent, bitangent, n).
    template <typename T>
    struct Basis {
        Vec3<T> tangent{};
        Vec3<T> bitangent{};
    };

    template <typename T>
    [[nodiscard]] constexpr Basis<T> orthonormalBasis(const Vec3<T>& n) noexcept
        requires std::is_floating_point_v<T>
    {
        CAD_ASSERT(n.isNormalized(T(1e-3)), "orthonormalBasis expects a normalized normal");

        const T sign = n.z >= T(0) ? T(1) : T(-1);
        const T a = T(-1) / (sign + n.z);
        const T b = n.x * n.y * a;

        return Basis<T>{
            Vec3<T>(T(1) + sign * n.x * n.x * a, sign* b, -sign * n.x),
                Vec3<T>(b, sign + n.y * n.y * a, -n.y)
        };
    }

    // -- stream output -------------------------------------------------------------------
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Vec3<T>& v) {
        return os << std::format("({}, {}, {})", v.x, v.y, v.z);
    }

} // namespace cad

// -- std::formatter / std::hash specializations --------------------------------
namespace std {

    template <typename T>
    struct formatter<cad::Vec3<T>, char> {
        cad::detail::FormatSpec spec;

        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }

        auto format(const cad::Vec3<T>& value, std::format_context& ctx) const {
            return cad::detail::emitBracketed(ctx, "(", ")", spec.view(),
                                              value.x, value.y, value.z);
        }
    };

    template <typename T>
    struct hash<cad::Vec3<T>> {
        [[nodiscard]] std::size_t operator()(const cad::Vec3<T>& v) const noexcept {
            std::size_t seed = std::hash<T>{}(v.x);
            seed = cad::hashCombine(seed, std::hash<T>{}(v.y));
            return cad::hashCombine(seed, std::hash<T>{}(v.z));
        }
    };

} // namespace std