#pragma once

#include <array>
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
#include "../support/Assertions.hpp"

namespace cad {

    template <typename T>
    struct Mat3 {
        static_assert(std::is_arithmetic_v<T>, "Mat3<T> requires an arithmetic T");

        using value_type = T;
        using Column = Vec3<T>;
        std::array<Vec3<T>, 3> cols{
            Vec3<T>(T(1), T(0), T(0)),
            Vec3<T>(T(0), T(1), T(0)),
            Vec3<T>(T(0), T(0), T(1)),
        };

        constexpr Mat3() noexcept = default;

        constexpr Mat3(const Vec3<T>& c0, const Vec3<T>& c1, const Vec3<T>& c2) noexcept
            : cols{c0, c1, c2} {
        }

        static constexpr Mat3 Identity() noexcept { return Mat3{}; }

        static constexpr Mat3 Zero() noexcept {
            return Mat3(Vec3<T>::Zero(), Vec3<T>::Zero(), Vec3<T>::Zero());
        }

        static constexpr Mat3 fromColumns(const Vec3<T>& c0, const Vec3<T>& c1,
                                          const Vec3<T>& c2) noexcept {
            return Mat3(c0, c1, c2);
        }

        // Rows are supplied in reading order: row0, row1, row2.
        static constexpr Mat3 fromRows(const Vec3<T>& r0, const Vec3<T>& r1,
                                       const Vec3<T>& r2) noexcept {
            return Mat3(Vec3<T>(r0.x, r1.x, r2.x),
                        Vec3<T>(r0.y, r1.y, r2.y),
                        Vec3<T>(r0.z, r1.z, r2.z));
        }

        static constexpr Mat3 fromDiagonal(T x, T y, T z) noexcept {
            return Mat3(Vec3<T>(x, T(0), T(0)),
                        Vec3<T>(T(0), y, T(0)),
                        Vec3<T>(T(0), T(0), z));
        }

        static constexpr Mat3 fromDiagonal(const Vec3<T>& d) noexcept {
            return fromDiagonal(d.x, d.y, d.z);
        }

        static constexpr Mat3 fromScale(T sx, T sy, T sz) noexcept {
            return fromDiagonal(sx, sy, sz);
        }

        static constexpr Mat3 fromScale(const Vec3<T>& s) noexcept {
            return fromDiagonal(s);
        }

        static constexpr Mat3 fromScale(T uniform) noexcept {
            return fromDiagonal(uniform, uniform, uniform);
        }

        // The columns ARE the images of the basis vectors.
        static constexpr Mat3 fromBasis(const Vec3<T>& xAxis, const Vec3<T>& yAxis,
                                        const Vec3<T>& zAxis) noexcept {
            return Mat3(xAxis, yAxis, zAxis);
        }

        [[nodiscard]] static Mat3 rotationX(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            const T c = std::cos(angleRadians);
            const T s = std::sin(angleRadians);
            return Mat3(Vec3<T>(T(1), T(0), T(0)),
                        Vec3<T>(T(0), c, s),
                        Vec3<T>(T(0), -s, c));
        }

        [[nodiscard]] static Mat3 rotationY(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            const T c = std::cos(angleRadians);
            const T s = std::sin(angleRadians);
            return Mat3(Vec3<T>(c, T(0), -s),
                        Vec3<T>(T(0), T(1), T(0)),
                        Vec3<T>(s, T(0), c));
        }

        [[nodiscard]] static Mat3 rotationZ(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            const T c = std::cos(angleRadians);
            const T s = std::sin(angleRadians);
            return Mat3(Vec3<T>(c, s, T(0)),
                        Vec3<T>(-s, c, T(0)),
                        Vec3<T>(T(0), T(0), T(1)));
        }

        // Rodrigues' rotation formula. axis must be unit length.
        [[nodiscard]] static Mat3 fromAxisAngle(const Vec3<T>& axis, T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(axis.isNormalized(T(1e-3)),
                       "Mat3::fromAxisAngle expects a normalized axis");

            const T c = std::cos(angleRadians);
            const T s = std::sin(angleRadians);
            const T k = T(1) - c;

            const T xx = axis.x * axis.x, xy = axis.x * axis.y, xz = axis.x * axis.z;
            const T yy = axis.y * axis.y, yz = axis.y * axis.z, zz = axis.z * axis.z;

            const T sx = axis.x * s, sy = axis.y * s, sz = axis.z * s;

            return Mat3(
                Vec3<T>(k * xx + c, k * xy + sz, k * xz - sy),
                Vec3<T>(k * xy - sz, k * yy + c, k * yz + sx),
                Vec3<T>(k * xz + sy, k * yz - sx, k * zz + c));
        }

        // element access
        [[nodiscard]] constexpr T& operator()(int col, int row) noexcept {
            CAD_ASSERT(col >= 0 && col < 3 && row >= 0 && row < 3,
                       "Mat3 index out of range");
            return cols[static_cast<std::size_t>(col)][row];
        }

        [[nodiscard]] constexpr const T& operator()(int col, int row) const noexcept {
            CAD_ASSERT(col >= 0 && col < 3 && row >= 0 && row < 3,
                       "Mat3 index out of range");
            return cols[static_cast<std::size_t>(col)][row];
        }

        [[nodiscard]] constexpr Vec3<T> getColumn(int col) const noexcept {
            CAD_ASSERT(col >= 0 && col < 3, "Mat3 column index out of range");
            return cols[static_cast<std::size_t>(col)];
        }

        [[nodiscard]] constexpr Vec3<T> getRow(int row) const noexcept {
            CAD_ASSERT(row >= 0 && row < 3, "Mat3 row index out of range");
            return Vec3<T>(cols[0][row], cols[1][row], cols[2][row]);
        }

        constexpr Mat3& operator*=(T s) noexcept {
            for(auto& c : cols) c *= s;
            return *this;
        }

        constexpr Mat3& operator/=(T s) noexcept {
            CAD_ASSERT(s != T(0), "Mat3 division by zero");
            for(auto& c : cols) c /= s;
            return *this;
        }

        [[nodiscard]] constexpr Mat3 transposed() const noexcept {
            return Mat3(getRow(0), getRow(1), getRow(2));
        }

        [[nodiscard]] constexpr T trace() const noexcept {
            return cols[0][0] + cols[1][1] + cols[2][2];
        }

        [[nodiscard]] constexpr Vec3<T> diagonal() const noexcept {
            return Vec3<T>(cols[0][0], cols[1][1], cols[2][2]);
        }

        [[nodiscard]] constexpr T minorOf(int skipRow, int skipCol) const noexcept {
            T m[4]{};
            int k = 0;
            for(int c = 0; c < 3; ++c) {
                if(c == skipCol) continue;
                for(int r = 0; r < 3; ++r) {
                    if(r == skipRow) continue;
                    m[k++] = cols[static_cast<std::size_t>(c)][r];
                }
            }
            return m[0] * m[3] - m[1] * m[2];
        }

        [[nodiscard]] constexpr T cofactor(int row, int col) const noexcept {
            const T m = minorOf(row, col);
            return ((row + col) & 1) == 0 ? m : -m;
        }

        [[nodiscard]] constexpr T determinant() const noexcept {
            return cols[0][0] * cofactor(0, 0)
                + cols[1][0] * cofactor(0, 1)
                + cols[2][0] * cofactor(0, 2);
        }

        [[nodiscard]] constexpr Mat3 inverted() const noexcept {
            const T det = determinant();
            CAD_VERIFY(det != T(0), "Mat3::inverted on singular matrix");

            Mat3 inv;
            // inverse(col c, row r) = cofactor(row c, col r) / det   (adjugate / det)
            for(int c = 0; c < 3; ++c) {
                for(int r = 0; r < 3; ++r) {
                    inv.cols[static_cast<std::size_t>(c)][r] = cofactor(c, r) / det;
                }
            }
            return inv;
        }

        [[nodiscard]] Mat3 orthonormalized() const noexcept
            requires std::is_floating_point_v<T>
        {
            Vec3<T> c0 = cols[0].normalizedOr(Vec3<T>::UnitX());

            Vec3<T> c1 = cols[1] - c0 * cols[1].dot(c0);
            c1 = c1.normalizedOr(Vec3<T>::UnitY() - c0 * c0.y); 

            Vec3<T> c2 = c0.cross(c1); 
            return Mat3(c0, c1, c2);
        }

        [[nodiscard]] constexpr bool isIdentity(
            T tol = default_abs_tol<T>) const noexcept
            requires std::is_floating_point_v<T>
        {
            return almostEqual(cols[0], Vec3<T>::UnitX(), tol, tol)
                && almostEqual(cols[1], Vec3<T>::UnitY(), tol, tol)
                && almostEqual(cols[2], Vec3<T>::UnitZ(), tol, tol);
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return cols[0].isFinite() && cols[1].isFinite() && cols[2].isFinite();
        }
    };

    using Mat3f = Mat3<float>;
    using Mat3d = Mat3<double>;

    // free operators 
    template <typename T>
    [[nodiscard]] constexpr Mat3<T> operator*(const Mat3<T>& a, const Mat3<T>& b) noexcept {
        Mat3<T> r;
        for(int i = 0; i < 3; ++i) {
            const Vec3<T>& bc = b.cols[static_cast<std::size_t>(i)];
            r.cols[static_cast<std::size_t>(i)] =
                a.cols[0] * bc.x + a.cols[1] * bc.y + a.cols[2] * bc.z;
        }
        return r;
    }

    template <typename T>
    [[nodiscard]] constexpr Vec3<T> operator*(const Mat3<T>& m, const Vec3<T>& v) noexcept {
        return m.cols[0] * v.x + m.cols[1] * v.y + m.cols[2] * v.z;
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> operator*(const Mat3<T>& m, T s) noexcept {
        return Mat3<T>(m.cols[0] * s, m.cols[1] * s, m.cols[2] * s);
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> operator*(T s, const Mat3<T>& m) noexcept {
        return m * s;
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> operator/(const Mat3<T>& m, T s) noexcept {
        CAD_ASSERT(s != T(0), "Mat3 division by zero");
        return Mat3<T>(m.cols[0] / s, m.cols[1] / s, m.cols[2] / s);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Mat3<T>& a, const Mat3<T>& b) noexcept {
        return a.cols == b.cols;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Mat3<T>& a, const Mat3<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Mat3<T>& a, const Mat3<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.cols[0], b.cols[0], relTol, absTol)
            && almostEqual(a.cols[1], b.cols[1], relTol, absTol)
            && almostEqual(a.cols[2], b.cols[2], relTol, absTol);
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> transpose(const Mat3<T>& m) noexcept {
        return m.transposed();
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> inverse(const Mat3<T>& m) noexcept {
        return m.inverted();
    }

    template <typename T>
    [[nodiscard]] constexpr T determinant(const Mat3<T>& m) noexcept {
        return m.determinant();
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Mat3<T>& m) {
        return os << std::format("{}", m);
    }

}

namespace std{
    template <typename T>
    struct formatter<cad::Mat3<T>, char> {
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Mat3<T>& m, std::format_context& ctx) const {
            std::string piece;
            piece.reserve(spec.view().size() + 3);
            piece += "{:";
            piece += spec.view();
            piece += '}';
            auto out = ctx.out();
            out = std::format_to(out, "[");
            for(int r = 0; r < 3; ++r) {
                out = std::format_to(out, "{}[", r == 0 ? "" : "\n ");
                for(int c = 0; c < 3; ++c) {
                    if(c != 0) out = std::format_to(out, ", ");
                    out = std::vformat_to(out, piece,
                                          std::make_format_args(m.cols[static_cast<std::size_t>(c)][r]));
                }
                out = std::format_to(out, "]");
            }
            return std::format_to(out, "]");
        }
    };
} 