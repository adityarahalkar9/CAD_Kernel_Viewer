
#pragma once

#include <cmath>
#include <cstddef>
#include <format>
#include <ostream>
#include <string>
#include <string_view>
#include <type_traits>

#include "Constants.hpp"
#include "Utils.hpp"
#include "Vec2.hpp"
#include "Vec3.hpp"
#include "Mat3.hpp"
#include "Mat4.hpp"
#include "../support/Assertions.hpp"

namespace cad {
    template <typename T>
    struct Quat {
        static_assert(std::is_floating_point_v<T>,
                      "Quat<T> requires a floating point T");

        using value_type = T;

        T w{1};
        T x{};
        T y{};
        T z{};

        constexpr Quat() noexcept = default;

        constexpr Quat(T iw, T ix, T iy, T iz) noexcept
            : w(iw), x(ix), y(iy), z(iz) {
        }

        static constexpr Quat Identity() noexcept { return Quat{}; }

        [[nodiscard]] static Quat fromAxisAngle(const Vec3<T>& axis, T angleRadians) noexcept {
            CAD_ASSERT(axis.isNormalized(T(1e-3)),
                       "Quat::fromAxisAngle expects a normalized axis");
            const T h = angleRadians / T(2);
            const T s = std::sin(h);
            return Quat(std::cos(h), axis.x * s, axis.y * s, axis.z * s);
        }

        [[nodiscard]] static Quat fromTwoVectors(const Vec3<T>& a, const Vec3<T>& b) noexcept {
            CAD_ASSERT(a.isNormalized(T(1e-3)) && b.isNormalized(T(1e-3)),
                       "Quat::fromTwoVectors expects normalized inputs");

            const T d = a.dot(b);

            if(d >= T(1) - T(64) * machine_eps<T>) {
                return Identity();
            }

            if(d <= T(-1) + T(64) * machine_eps<T>) {
                const Basis<T> basis = orthonormalBasis(a);
                return fromAxisAngle(basis.tangent, pi_v<T>);
            }

            const Vec3<T> c = a.cross(b);
            Quat q(T(1) + d, c.x, c.y, c.z);
            return q.normalized();
        }

        [[nodiscard]] static Quat fromEulerXYZ(T xRad, T yRad, T zRad) noexcept {
            const Quat qx = fromAxisAngle(Vec3<T>::UnitX(), xRad);
            const Quat qy = fromAxisAngle(Vec3<T>::UnitY(), yRad);
            const Quat qz = fromAxisAngle(Vec3<T>::UnitZ(), zRad);
            return qz * qy * qx;
        }

        [[nodiscard]] static Quat fromMat3(const Mat3<T>& m) noexcept {
            const T trace = m.cols[0][0] + m.cols[1][1] + m.cols[2][2];

            if(trace > T(0)) {
                const T s = cad::sqrt(trace + T(1)) * T(2); // s = 4w
                return Quat(s / T(4),
                            (m.cols[2][1] - m.cols[1][2]) / s,
                            (m.cols[0][2] - m.cols[2][0]) / s,
                            (m.cols[1][0] - m.cols[0][1]) / s);
            }
            if(m.cols[0][0] > m.cols[1][1] && m.cols[0][0] > m.cols[2][2]) {
                const T s = cad::sqrt(T(1) + m.cols[0][0] - m.cols[1][1] - m.cols[2][2]) * T(2);
                return Quat((m.cols[2][1] - m.cols[1][2]) / s,
                            s / T(4),
                            (m.cols[0][1] + m.cols[1][0]) / s,
                            (m.cols[0][2] + m.cols[2][0]) / s);
            }
            if(m.cols[1][1] > m.cols[2][2]) {
                const T s = cad::sqrt(T(1) + m.cols[1][1] - m.cols[0][0] - m.cols[2][2]) * T(2);
                return Quat((m.cols[0][2] - m.cols[2][0]) / s,
                            (m.cols[0][1] + m.cols[1][0]) / s,
                            s / T(4),
                            (m.cols[1][2] + m.cols[2][1]) / s);
            }
            const T s = cad::sqrt(T(1) + m.cols[2][2] - m.cols[0][0] - m.cols[1][1]) * T(2);
            return Quat((m.cols[1][0] - m.cols[0][1]) / s,
                        (m.cols[0][2] + m.cols[2][0]) / s,
                        (m.cols[1][2] + m.cols[2][1]) / s,
                        s / T(4));
        }

        [[nodiscard]] static Quat fromMat4(const Mat4<T>& m) noexcept {
            return fromMat3(m.upperLeft());
        }

        // basic ops
        [[nodiscard]] constexpr T dot(const Quat& r) const noexcept {
            return w * r.w + x * r.x + y * r.y + z * r.z;
        }

        [[nodiscard]] constexpr T normSquared() const noexcept {
            return w * w + x * x + y * y + z * z;
        }

        [[nodiscard]] T norm() const noexcept {
            return cad::sqrt(normSquared());
        }

        [[nodiscard]] constexpr bool isZero() const noexcept {
            return w == T(0) && x == T(0) && y == T(0) && z == T(0);
        }

        [[nodiscard]] constexpr Quat conjugate() const noexcept {
            return Quat(w, -x, -y, -z);
        }

        [[nodiscard]] constexpr Quat inverted() const noexcept {
            CAD_ASSERT(normSquared() > T(0), "Quat::inverted on zero quaternion");
            const T invNorm2 = T(1) / normSquared();
            return Quat(w * invNorm2, -x * invNorm2, -y * invNorm2, -z * invNorm2);
        }

        [[nodiscard]] constexpr Quat normalized() const noexcept {
            CAD_ASSERT(normSquared() > T(0), "Quat::normalized on zero quaternion");
            const T inv = T(1) / cad::sqrt(normSquared());
            return Quat(w * inv, x * inv, y * inv, z * inv);
        }

        [[nodiscard]] constexpr bool isNormalized(
            T tol = GeometryTolerance<T>::unit_deviation) const noexcept {
            return almostZero(normSquared() - T(1), tol);
        }

        // rotation extraction 
        [[nodiscard]] T angle() const noexcept {
            Quat q = w < T(0) ? -*this : *this;
            const T cw = q.w > T(1) ? T(1) : (q.w < T(-1) ? T(-1) : q.w);
            return T(2) * std::acos(cw);
        }

        // Rotation axis
        [[nodiscard]] Vec3<T> axis() const noexcept {
            Quat q = w < T(0) ? -*this : *this;               // shortest form
            const Vec3<T> v(q.x, q.y, q.z);
            const T s2 = T(1) - q.w * q.w;
            if(s2 <= T(64) * machine_eps<T> *T(64) * machine_eps<T>) {
                return Vec3<T>::UnitX();
            }
            return v / cad::sqrt(s2);
        }

        // rotating vectors 
        [[nodiscard]] constexpr Vec3<T> rotate(const Vec3<T>& p) const noexcept {
            CAD_ASSERT(isNormalized(T(1e-3)), "Quat::rotate expects a unit quaternion");
            const Vec3<T> v(x, y, z);
            const Vec3<T> t = v.cross(p) * T(2);
            return p + t * w + v.cross(t);
        }

        // matrix conversion
        [[nodiscard]] constexpr Mat3<T> toMat3() const noexcept {
            CAD_ASSERT(isNormalized(T(1e-3)), "Quat::toMat3 expects a unit quaternion");

            const T xx = x * x, yy = y * y, zz = z * z;
            const T xy = x * y, xz = x * z, yz = y * z;
            const T wx = w * x, wy = w * y, wz = w * z;

            return Mat3<T>(
                Vec3<T>(T(1) - T(2) * (yy + zz), T(2) * (xy + wz), T(2) * (xz - wy)),
                Vec3<T>(T(2) * (xy - wz), T(1) - T(2) * (xx + zz), T(2) * (yz + wx)),
                Vec3<T>(T(2) * (xz + wy), T(2) * (yz - wx), T(1) - T(2) * (xx + yy)));
        }

        [[nodiscard]] constexpr Mat4<T> toMat4() const noexcept {
            return Mat4<T>::fromMat3(toMat3(), Vec3<T>::Zero());
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return std::isfinite(w) && std::isfinite(x)
                && std::isfinite(y) && std::isfinite(z);
        }
    };

    // aliases
    using Quatf = Quat<float>;
    using Quatd = Quat<double>;

    // Hamilton product: (a * b) 
    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator*(const Quat<T>& a, const Quat<T>& b) noexcept {
        return Quat<T>(
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w);
    }

    template <typename T>
    constexpr Quat<T>& operator*=(Quat<T>& a, const Quat<T>& b) noexcept {
        a = a * b;
        return a;
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator-(const Quat<T>& q) noexcept {
        return Quat<T>(-q.w, -q.x, -q.y, -q.z);
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator+(const Quat<T>& a, const Quat<T>& b) noexcept {
        return Quat<T>(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z);
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator-(const Quat<T>& a, const Quat<T>& b) noexcept {
        return Quat<T>(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z);
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator*(const Quat<T>& q, T s) noexcept {
        return Quat<T>(q.w * s, q.x * s, q.y * s, q.z * s);
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator*(T s, const Quat<T>& q) noexcept {
        return q * s;
    }

    template <typename T>
    [[nodiscard]] constexpr Quat<T> operator/(const Quat<T>& q, T s) noexcept {
        CAD_ASSERT(s != T(0), "Quat division by zero");
        return Quat<T>(q.w / s, q.x / s, q.y / s, q.z / s);
    }

    // comparison
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Quat<T>& a, const Quat<T>& b) noexcept {
        return a.w == b.w && a.x == b.x && a.y == b.y && a.z == b.z;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Quat<T>& a, const Quat<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Quat<T>& a, const Quat<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        return almostEqual(a.w, b.w, relTol, absTol)
            && almostEqual(a.x, b.x, relTol, absTol)
            && almostEqual(a.y, b.y, relTol, absTol)
            && almostEqual(a.z, b.z, relTol, absTol);
    }

    template <typename T>
    [[nodiscard]] constexpr T dot(const Quat<T>& a, const Quat<T>& b) noexcept {
        return a.dot(b);
    }

    // interpolation 
    template <typename T>
    [[nodiscard]] Quat<T> nlerp(const Quat<T>& a, const Quat<T>& b, T t) noexcept {
        CAD_ASSERT(t >= T(0) && t <= T(1), "nlerp: t outside [0, 1]");
        Quat<T> bb = b;
        if(a.dot(bb) < T(0)) {
            bb = -bb; // shortest arc
        }
        return (a + (bb - a) * t).normalized();
    }

    template <typename T>
    [[nodiscard]] Quat<T> slerp(const Quat<T>& a, const Quat<T>& b, T t) noexcept {
        CAD_ASSERT(t >= T(0) && t <= T(1), "slerp: t outside [0, 1]");

        T d = a.dot(b);
        Quat<T> bb = b;
        if(d < T(0)) {
            bb = -bb;
            d = -d;
        }

        constexpr T kNearParallel = T(0.9995);
        if(d >= kNearParallel) {
            return nlerp(a, bb, t);
        }

        const T theta = std::acos(d);
        const T sinTheta = std::sin(theta);

        const T w0 = std::sin((T(1) - t) * theta) / sinTheta;
        const T w1 = std::sin(t * theta) / sinTheta;

        return (a * w0) + (bb * w1);
    }

    // stream output
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Quat<T>& q) {
        return os << std::format("{}", q);
    }

}

// std::formatter 
namespace std {

    template <typename T>
    struct formatter<cad::Quat<T>, char> {
        cad::detail::FormatSpec spec;

        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }

        // (w, x, y, z) - w first, matching storage and constructor order.
        auto format(const cad::Quat<T>& q, std::format_context& ctx) const {
            return cad::detail::emitBracketed(ctx, "(", ")", spec.view(),
                                              q.w, q.x, q.y, q.z);
        }
    };
} 