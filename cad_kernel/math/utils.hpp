// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : math/Utils.hpp
//  Purpose : Numeric utilities used across the kernel
//
//  Notable: cad::sqrt is usable BOTH at compile time (Newton iteration with
//  guaranteed termination) and at runtime (std::sqrt). C++20 has no constexpr
//  std::sqrt; this removes that ceiling so normalize()/length() are genuine
//  constexpr - verified again by unit tests via static_assert.
// ============================================================================
#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "Constants.hpp"
#include "../support/Assertions.hpp"

namespace cad {

    // ---------------------------------------------------------------------------
    // constexpr sqrt
    // ---------------------------------------------------------------------------
    namespace detail {

        // 2^k for integer k (std::pow / std::ldexp are not constexpr in C++20).
        template <typename T>
        [[nodiscard]] constexpr T pow2Int(int k) {
            T m = T(1);
            for(int i = 0; i < k; ++i)   m *= T(2);
            for(int i = 0; i > k; --i)   m *= T(0.5);
            return m;
        }

        // Newton-Raphson sqrt with guaranteed termination (compile-time path only).
        //  * argument scaled into [1, 4) so the iteration starts close to the root
        //    (quadratic convergence: <= ~8 iterations for double)
        //  * after the first step the sequence decreases monotonically, so
        //    "next >= current" is a safe fixed-point / rounding-oscillation stop
        template <typename T>
        [[nodiscard]] constexpr T sqrtNewton(T x) {
            if(x != x)                    return std::numeric_limits<T>::quiet_NaN();
            if(x == T(0))                 return T(0);          // covers -0.0 too
            if(x < T(0))                 return std::numeric_limits<T>::quiet_NaN();
            if(x == cad::infinity<T>())   return cad::infinity<T>();

            // sqrt(x) = sqrt(s) * 2^k  with s in [1, 4)
            T   s = x;
            int k = 0;
            while(s >= T(4)) { s *= T(0.25); ++k; }
            while(s < T(1)) { s *= T(4);   --k; }

            T r = s;   // in [1,4) -> sqrt(s) in [1,2)
            for(int i = 0; i < 64; ++i) {
                const T next = (r + s / r) / T(2);
                if(next >= r) break;
                r = next;
            }
            return r * pow2Int<T>(k);
        }

    } // namespace detail

    // sqrt for compile time AND runtime.
    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T sqrt(T x) {
        if(std::is_constant_evaluated()) {
            return detail::sqrtNewton(x);
        }
        return std::sqrt(x);
    }

    // ---------------------------------------------------------------------------
    // Small numeric helpers
    // ---------------------------------------------------------------------------
    template <typename T>
    [[nodiscard]] constexpr T sqr(T x) noexcept { return x * x; }

    // -1 / 0 / +1 (never NaN-propagating; sign(-0.0) == 0).
    template <typename T>
    [[nodiscard]] constexpr T sign(T v) noexcept {
        return T(T(0) < v) - T(v < T(0));
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T lerp(T a, T b, T t) noexcept {
        // t is intentionally unconstrained: t outside [0,1] extrapolates.
        return a + (b - a) * t;
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T saturate(T v) noexcept {
        return v < T(0) ? T(0) : (v > T(1) ? T(1) : v);
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T smoothstep(T edge0, T edge1, T x) noexcept {
        CAD_ASSERT(edge0 != edge1, "smoothstep: edge0 == edge1");
        const T t = saturate((x - edge0) / (edge1 - edge0));
        return t * t * (T(3) - T(2) * t);
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T radians(T degrees) noexcept {
        return degrees * deg_to_rad_v<T>;
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr T degrees(T radians) noexcept {
        return radians * rad_to_deg_v<T>;
    }

    // ---------------------------------------------------------------------------
    // Approximate equality (floats)
    //
    //  * NaN never compares equal to anything (including itself)
    //  * +/-inf compare equal only to themselves (handled by the a == b fast path)
    //  * machine-scaled defaults: 64 eps ~ 1.4e-14 for double
    // ---------------------------------------------------------------------------
    template <typename T>
        requires std::is_floating_point_v<T>
    inline constexpr T default_rel_tol = T(64) * machine_eps<T>;

    template <typename T>
        requires std::is_floating_point_v<T>
    inline constexpr T default_abs_tol = T(64) * machine_eps<T>;

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr bool almostEqual(T a, T b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        if(a == b) {
            return true;
        }
        const T d = a > b ? a - b : b - a;
        if(d <= absTol) {
            return true;
        }
        const T aa = a < T(0) ? -a : a;
        const T ab = b < T(0) ? -b : b;
        const T m = aa > ab ? aa : ab;
        return d <= relTol * m;
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] constexpr bool almostZero(T v, T absTol = default_abs_tol<T>) noexcept {
        return (v < T(0) ? -v : v) <= absTol;
    }

    // ---------------------------------------------------------------------------
    // Hashing (std::hash glue for the vector types, and mesh weld keys later)
    // ---------------------------------------------------------------------------
    [[nodiscard]] constexpr std::size_t hashCombine(std::size_t seed, std::size_t value) noexcept {
        return seed ^ (value + 0x9e3779b97f4a7c15ull
                       + (seed << 6) + (seed >> 2));
    }

    // ---------------------------------------------------------------------------
    // Quantization (runtime - used by vertex welding in mesh/MeshAlgorithms)
    //
    // Values are snapped to the CENTER of a grid cell. Two points closer than
    // `grid` can still land in adjacent cells, so the welder probes neighbouring
    // cells (27 buckets in 3D) rather than trusting a single bucket lookup.
    // ---------------------------------------------------------------------------
    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] inline T quantize(T v, T grid) noexcept {
        CAD_ASSERT(grid > T(0), "quantize: grid must be > 0");
        return std::floor(v / grid + T(0.5)) * grid;
    }

    template <typename T>
        requires std::is_floating_point_v<T>
    [[nodiscard]] inline std::int64_t bucketIndex(T v, T grid) noexcept {
        CAD_ASSERT(grid > T(0), "bucketIndex: grid must be > 0");
        // Assumes |v / grid| < 9.2e18 (always true for sane CAD coordinates).
        return static_cast<std::int64_t>(std::floor(v / grid + T(0.5)));
    }

} // namespace cad