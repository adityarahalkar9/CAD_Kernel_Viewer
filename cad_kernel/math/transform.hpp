
#pragma once

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
#include "Vec4.hpp"
#include "Mat3.hpp"
#include "Mat4.hpp"
#include "Quaternion.hpp"
#include "../support/Assertions.hpp"

namespace cad {

    template <typename T>
    struct Transform {
        static_assert(std::is_floating_point_v<T>,
                      "Transform<T> requires a floating point T");

        using value_type = T;
        Vec3<T> position{};
        Quat<T> rotation{};
        Vec3<T> scale{T(1), T(1), T(1)};

        constexpr Transform() noexcept = default;

        constexpr Transform(const Vec3<T>& position, const Quat<T>& rotation,
                            const Vec3<T>& scale) noexcept
            : position(position), rotation(rotation), scale(scale) {
        }

        static constexpr Transform Identity() noexcept { return Transform{}; }

        static constexpr Transform fromTranslation(const Vec3<T>& t) noexcept {
            return Transform(t, Quat<T>::Identity(), Vec3<T>::One());
        }

        static constexpr Transform fromRotation(const Quat<T>& q) noexcept {
            return Transform(Vec3<T>::Zero(), q, Vec3<T>::One());
        }

        static constexpr Transform fromScale(const Vec3<T>& s) noexcept {
            return Transform(Vec3<T>::Zero(), Quat<T>::Identity(), s);
        }

        static constexpr Transform fromScale(T uniform) noexcept {
            return fromScale(Vec3<T>(uniform, uniform, uniform));
        }

        // M = T * R * S.
        [[nodiscard]] constexpr Mat4<T> toMatrix() const noexcept {
            Mat3<T> lin = rotation.toMat3();
            lin.cols[0] = lin.cols[0] * scale.x;
            lin.cols[1] = lin.cols[1] * scale.y;
            lin.cols[2] = lin.cols[2] * scale.z;
            return Mat4<T>::fromMat3(lin, position);
        }

        [[nodiscard]] static Transform fromMatrix(const Mat4<T>& m) noexcept {
            Mat3<T> lin = m.upperLeft();

            if(lin.determinant() < T(0)) {
                lin.cols[0] = -lin.cols[0]; // fold mirror into x
            }

            const Vec3<T> sx = lin.cols[0].lengthSquared() > T(0)
                ? Vec3<T>(cad::sqrt(lin.cols[0].lengthSquared()), T(0), T(0))
                : Vec3<T>::Zero();
            const T sy = lin.cols[1].length();
            const T sz = lin.cols[2].length();

            CAD_VERIFY(sx.x > T(0) && sy > T(0) && sz > T(0),
                       "Transform::fromMatrix: degenerate (zero) scale axis");

            const Vec3<T> s(sx.x, sy, sz);

            Mat3<T> rot = lin;
            rot.cols[0] = rot.cols[0] / s.x;
            rot.cols[1] = rot.cols[1] / s.y;
            rot.cols[2] = rot.cols[2] / s.z;

            return Transform(m.getTranslation(), Quat<T>::fromMat3(rot), s);
        }

        // direct point math 
        [[nodiscard]] constexpr Vec3<T> transformPoint(const Vec3<T>& p) const noexcept {
            return position + rotation.rotate(cwMul(scale, p));
        }

        // Scale + rotate, no translation.
        [[nodiscard]] constexpr Vec3<T> transformDirection(const Vec3<T>& d) const noexcept {
            return rotation.rotate(cwMul(scale, d));
        }

        // inverse
        [[nodiscard]] constexpr Transform inverted() const noexcept {
            CAD_ASSERT(hasUniformScale(), "Transform::inverted requires uniform scale");

            const Quat<T> invRot = rotation.conjugate(); // unit quaternion
            const T       invS = T(1) / scale.x;

            Transform r;
            r.rotation = invRot;
            r.scale = Vec3<T>(invS, invS, invS);
            r.position = invRot.rotate(position) * (-invS);
            return r;
        }

        // General exact inverse for any TRS (goes through the matrix).
        [[nodiscard]] constexpr Mat4<T> toInverseMatrix() const noexcept {
            return toMatrix().inverted();
        }

        // composition 
        [[nodiscard]] static Transform compose(const Transform& parent,
                                               const Transform& child) noexcept {
            return fromMatrix(parent.toMatrix() * child.toMatrix());
        }

        // queries 
        [[nodiscard]] constexpr bool hasUniformScale(
            T tol = default_rel_tol<T>) const noexcept {
            return almostEqual(scale.x, scale.y, tol, tol)
                && almostEqual(scale.y, scale.z, tol, tol);
        }

        // Any non-positive scale component collapses or mirrors geometry.
        [[nodiscard]] constexpr bool hasNonPositiveScale() const noexcept {
            return scale.x <= T(0) || scale.y <= T(0) || scale.z <= T(0);
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return position.isFinite() && scale.isFinite() && rotation.isFinite();
        }
    };

    // aliases 
    using Transformf = Transform<float>;
    using Transformd = Transform<double>;

    // comparison 
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Transform<T>& a,
                                            const Transform<T>& b) noexcept {
        return a.position == b.position && a.rotation == b.rotation && a.scale == b.scale;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Transform<T>& a,
                                            const Transform<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Transform<T>& a, const Transform<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        return almostEqual(a.position, b.position, relTol, absTol)
            && almostEqual(a.rotation, b.rotation, relTol, absTol)
            && almostEqual(a.scale, b.scale, relTol, absTol);
    }

    // interpolation 
    template <typename T>
    [[nodiscard]] Transform<T> lerp(const Transform<T>& a, const Transform<T>& b, T t) noexcept {
        return Transform<T>(cad::lerp(a.position, b.position, t),
                            slerp(a.rotation, b.rotation, t),
                            cad::lerp(a.scale, b.scale, t));
    }

    // stream output
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Transform<T>& t) {
        return os << std::format("{}", t);
    }

} 

// std::formatter 
namespace std{
    template <typename T>
    struct formatter<cad::Transform<T>, char> {
        cad::detail::FormatSpec spec;

        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }

        auto format(const cad::Transform<T>& value, std::format_context& ctx) const {
            // Reuses the Vec3/Quat formatters with the same spec.
            std::string piece;
            piece.reserve(spec.view().size() + 3);
            piece += "{:";
            piece += spec.view();
            piece += '}';

            const std::string text = std::vformat(
                "[p " + piece + ", q " + piece + ", s " + piece + "]",
                std::make_format_args(value.position, value.rotation, value.scale));

            return std::format_to(ctx.out(), "{}", text);
        }
    };
} 