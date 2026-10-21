// ============================================================================
//  cad_kernel - header-only CAD geometry kernel
//  File    : math/Constants.hpp
//  Purpose : Mathematical constants and kernel-wide tolerance defaults
//
//  Variable templates give per-type constants; the non-template double
//  constants are the everyday spelling ground (pi, twoPi, ...).
//  GeometryTolerance carries the documented "magic numbers" the kernel uses
//  when callers do not pass explicit tolerances.
// ============================================================================
#pragma once

#include <limits>
#include <type_traits>

namespace cad {

    // ---------------------------------------------------------------------------
    // Constants (variable templates + double spellings)
    // ---------------------------------------------------------------------------
    template <typename T>
    inline constexpr T pi_v = T(3.141592653589793238462643383279502884L);

    template <typename T> inline constexpr T two_pi_v = pi_v<T> *T(2);
    template <typename T> inline constexpr T half_pi_v = pi_v<T> / T(2);
    template <typename T> inline constexpr T inv_pi_v = T(1) / pi_v<T>;
    template <typename T> inline constexpr T inv_two_pi_v = T(1) / two_pi_v<T>;
    template <typename T> inline constexpr T deg_to_rad_v = pi_v<T> / T(180);
    template <typename T> inline constexpr T rad_to_deg_v = T(180) / pi_v<T>;
    template <typename T> inline constexpr T e_v = T(2.718281828459045235360287471352662498L);
    template <typename T> inline constexpr T sqrt2_v = T(1.414213562373095048801688724209698079L);

    inline constexpr double pi = pi_v<double>;
    inline constexpr double two_pi = two_pi_v<double>;
    inline constexpr double half_pi = half_pi_v<double>;
    inline constexpr double inv_pi = inv_pi_v<double>;
    inline constexpr double inv_two_pi = inv_two_pi_v<double>;
    inline constexpr double deg_to_rad = deg_to_rad_v<double>;
    inline constexpr double rad_to_deg = rad_to_deg_v<double>;
    inline constexpr double e = e_v<double>;
    inline constexpr double sqrt2 = sqrt2_v<double>;

    // Machine epsilon, spelled once.
    template <typename T>
    inline constexpr T machine_eps = std::numeric_limits<T>::epsilon();

    template <typename T>
        requires std::is_arithmetic_v<T>
    [[nodiscard]] constexpr T infinity() noexcept {
        return std::numeric_limits<T>::infinity();
    }

    // ---------------------------------------------------------------------------
    // Kernel tolerance defaults
    //
    // Tuned for the double-precision kernel (the float instantiations exist for
    // GPU-side convenience; pass explicit tolerances where it matters).
    // ---------------------------------------------------------------------------
    template <typename T>
    struct GeometryTolerance {
        static constexpr T weld = T(1e-6);  // vertex welding distance
        static constexpr T coincidence = T(1e-9);  // "same point" in tests
        static constexpr T unit_deviation = T(64) * machine_eps<T>; // |len - 1| check
        static constexpr T parallel = T(1e-8); // determinant ~ 0
        static constexpr T point_on_plane = T(1e-9);
    };

    static_assert(GeometryTolerance<double>::weld > 0.0);
    static_assert(GeometryTolerance<double>::unit_deviation > 0.0);

} // namespace cad