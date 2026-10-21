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
#include "Vec4.hpp"
#include "Mat3.hpp"
#include "../support/Assertions.hpp"

namespace cad {

    template <typename T>
    struct Mat4 {
        static_assert(std::is_arithmetic_v<T>, "Mat4<T> requires an arithmetic T");

        using value_type = T;
        using Column = Vec4<T>;

        // data
        std::array<Vec4<T>, 4> cols{
            Vec4<T>(T(1), T(0), T(0), T(0)),
            Vec4<T>(T(0), T(1), T(0), T(0)),
            Vec4<T>(T(0), T(0), T(1), T(0)),
            Vec4<T>(T(0), T(0), T(0), T(1)),
        };

        constexpr Mat4() noexcept = default;

        constexpr Mat4(const Vec4<T>& c0, const Vec4<T>& c1,
                       const Vec4<T>& c2, const Vec4<T>& c3) noexcept
            : cols{c0, c1, c2, c3} {
        }

        static constexpr Mat4 Identity() noexcept { return Mat4{}; }

        static constexpr Mat4 Zero() noexcept {
            return Mat4(Vec4<T>::Zero(), Vec4<T>::Zero(),
                        Vec4<T>::Zero(), Vec4<T>::Zero());
        }

        static constexpr Mat4 fromColumns(const Vec4<T>& c0, const Vec4<T>& c1,
                                          const Vec4<T>& c2, const Vec4<T>& c3) noexcept {
            return Mat4(c0, c1, c2, c3);
        }

        // affine building blocks 
        static constexpr Mat4 translation(const Vec3<T>& t) noexcept {
            return Mat4(Vec4<T>(T(1), T(0), T(0), T(0)),
                        Vec4<T>(T(0), T(1), T(0), T(0)),
                        Vec4<T>(T(0), T(0), T(1), T(0)),
                        Vec4<T>(t.x, t.y, t.z, T(1)));
        }

        static constexpr Mat4 translation(T x, T y, T z) noexcept {
            return translation(Vec3<T>(x, y, z));
        }

        static constexpr Mat4 scale(const Vec3<T>& s) noexcept {
            return Mat4(Vec4<T>(s.x, T(0), T(0), T(0)),
                        Vec4<T>(T(0), s.y, T(0), T(0)),
                        Vec4<T>(T(0), T(0), s.z, T(0)),
                        Vec4<T>(T(0), T(0), T(0), T(1)));
        }

        static constexpr Mat4 scale(T uniform) noexcept {
            return scale(Vec3<T>(uniform, uniform, uniform));
        }

        static constexpr Mat4 scale(T sx, T sy, T sz) noexcept {
            return scale(Vec3<T>(sx, sy, sz));
        }

        // Linear part + translation (bottom row 0,0,0,1). THE TRS assembler.
        static constexpr Mat4 fromMat3(const Mat3<T>& linear, const Vec3<T>& t) noexcept {
            return Mat4(
                Vec4<T>(linear.cols[0].x, linear.cols[0].y, linear.cols[0].z, T(0)),
                Vec4<T>(linear.cols[1].x, linear.cols[1].y, linear.cols[1].z, T(0)),
                Vec4<T>(linear.cols[2].x, linear.cols[2].y, linear.cols[2].z, T(0)),
                Vec4<T>(t.x, t.y, t.z, T(1)));
        }

        // Rotations reuse Mat3 (runtime only, cos/sin).
        [[nodiscard]] static Mat4 rotationX(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            return fromMat3(Mat3<T>::rotationX(angleRadians), Vec3<T>::Zero());
        }

        [[nodiscard]] static Mat4 rotationY(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            return fromMat3(Mat3<T>::rotationY(angleRadians), Vec3<T>::Zero());
        }

        [[nodiscard]] static Mat4 rotationZ(T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            return fromMat3(Mat3<T>::rotationZ(angleRadians), Vec3<T>::Zero());
        }

        [[nodiscard]] static Mat4 fromAxisAngle(const Vec3<T>& axis, T angleRadians) noexcept
            requires std::is_floating_point_v<T>
        {
            return fromMat3(Mat3<T>::fromAxisAngle(axis, angleRadians), Vec3<T>::Zero());
        }

        [[nodiscard]] static Mat4 perspective(T fovYRadians, T aspect,
                                              T near, T far) noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(near > T(0), "perspective: near must be > 0");
            CAD_ASSERT(far > near, "perspective: far must be > near");
            CAD_ASSERT(aspect > T(0), "perspective: aspect must be > 0");
            CAD_ASSERT(fovYRadians > T(0) && fovYRadians < pi_v<T>,
                       "perspective: fovY must be in (0, pi)");

            const T f = T(1) / std::tan(fovYRadians / T(2));

            return Mat4(
                Vec4<T>(f / aspect, T(0), T(0), T(0)),
                Vec4<T>(T(0), -f, T(0), T(0)),   // Vulkan y-flip
                Vec4<T>(T(0), T(0), far / (near - far), -T(1)),
                Vec4<T>(T(0), T(0), far * near / (near - far), T(0)));
        }

        [[nodiscard]] static Mat4 ortho(T left, T right, T bottom, T top,
                                        T near, T far) noexcept
            requires std::is_floating_point_v<T>
        {
            CAD_ASSERT(right > left, "ortho: right must be > left");
            CAD_ASSERT(top > bottom, "ortho: top must be > bottom");
            CAD_ASSERT(far > near && near >= T(0), "ortho: need far > near >= 0");

            return Mat4(
                Vec4<T>(T(2) / (right - left), T(0), T(0), T(0)),
                Vec4<T>(T(0), T(-2) / (top - bottom), T(0), T(0)),            // Vulkan y-flip
                Vec4<T>(T(0), T(0), T(-1) / (far - near), T(0)),
                Vec4<T>(-(left + right) / (right - left),
                        (top + bottom) / (top - bottom),
                        -near / (far - near),
                        T(1)));
        }

        [[nodiscard]] static Mat4 lookAt(const Vec3<T>& eye, const Vec3<T>& target,
                                         const Vec3<T>& up) noexcept
            requires std::is_floating_point_v<T>
        {
            const Vec3<T> forward = (target - eye).normalized();

            const Vec3<T> rightRaw = forward.cross(up);
            CAD_VERIFY(rightRaw.lengthSquared() > GeometryTolerance<T>::parallel *
                       up.lengthSquared() * forward.lengthSquared(),
                       "Mat4::lookAt: forward and up are (anti)parallel");
            const Vec3<T> right = rightRaw.normalized();

            const Vec3<T> upOrtho = right.cross(forward);

            return Mat4(
                Vec4<T>(right.x, upOrtho.x, -forward.x, T(0)),
                Vec4<T>(right.y, upOrtho.y, -forward.y, T(0)),
                Vec4<T>(right.z, upOrtho.z, -forward.z, T(0)),
                Vec4<T>(-right.dot(eye), -upOrtho.dot(eye), forward.dot(eye), T(1)));
        }

        // element access
        [[nodiscard]] constexpr T& operator()(int col, int row) noexcept {
            CAD_ASSERT(col >= 0 && col < 4 && row >= 0 && row < 4,
                       "Mat4 index out of range");
            return cols[static_cast<std::size_t>(col)][row];
        }

        [[nodiscard]] constexpr const T& operator()(int col, int row) const noexcept {
            CAD_ASSERT(col >= 0 && col < 4 && row >= 0 && row < 4,
                       "Mat4 index out of range");
            return cols[static_cast<std::size_t>(col)][row];
        }

        [[nodiscard]] constexpr Vec4<T> getColumn(int col) const noexcept {
            CAD_ASSERT(col >= 0 && col < 4, "Mat4 column index out of range");
            return cols[static_cast<std::size_t>(col)];
        }

        constexpr Mat4& setColumn(int col, const Vec4<T>& v) noexcept {
            CAD_ASSERT(col >= 0 && col < 4, "Mat4 column index out of range");
            cols[static_cast<std::size_t>(col)] = v;
            return *this;
        }

        [[nodiscard]] constexpr Vec4<T> getRow(int row) const noexcept {
            CAD_ASSERT(row >= 0 && row < 4, "Mat4 row index out of range");
            return Vec4<T>(cols[0][row], cols[1][row], cols[2][row], cols[3][row]);
        }

        // Linear (rotation/scale) 3x3 part.
        [[nodiscard]] constexpr Mat3<T> upperLeft() const noexcept {
            return Mat3<T>(Vec3<T>(cols[0].xyz()),
                           Vec3<T>(cols[1].xyz()),
                           Vec3<T>(cols[2].xyz()));
        }

        [[nodiscard]] constexpr Vec3<T> getTranslation() const noexcept {
            return cols[3].xyz();
        }

        // point / direction transforms
        [[nodiscard]] constexpr Vec3<T> transformPoint(const Vec3<T>& p) const noexcept {
            const Vec3<T> c0 = cols[0].xyz();
            const Vec3<T> c1 = cols[1].xyz();
            const Vec3<T> c2 = cols[2].xyz();
            const Vec3<T> c3 = cols[3].xyz();
            return c0 * p.x + c1 * p.y + c2 * p.z + c3;
        }

        // General path with homogeneous divide
        [[nodiscard]] constexpr Vec3<T> transformCoord(const Vec3<T>& p) const noexcept
            requires std::is_floating_point_v<T>
        {
            return ((*this) * Vec4<T>::point(p)).perspectiveDivide();
        }

        [[nodiscard]] constexpr Vec3<T> transformDir(const Vec3<T>& d) const noexcept {
            const Vec3<T> c0 = cols[0].xyz();
            const Vec3<T> c1 = cols[1].xyz();
            const Vec3<T> c2 = cols[2].xyz();
            return c0 * d.x + c1 * d.y + c2 * d.z;
        }

        // arithmetic 
        constexpr Mat4& operator*=(T s) noexcept {
            for(auto& c : cols) c *= s;
            return *this;
        }
        constexpr Mat4& operator/=(T s) noexcept {
            CAD_ASSERT(s != T(0), "Mat4 division by zero");
            for(auto& c : cols) c /= s;
            return *this;
        }
        [[nodiscard]] constexpr Mat4 transposed() const noexcept {
            return Mat4(getRow(0), getRow(1), getRow(2), getRow(3));
        }
        constexpr Mat4& dropTranslation() noexcept {
            cols[3] = Vec4<T>(T(0), T(0), T(0), T(1));
            return *this;
        }

        // queries 
        [[nodiscard]] constexpr T minorOf(int skipRow, int skipCol) const noexcept {
            T m[9]{};
            int k = 0;
            for(int c = 0; c < 4; ++c) {
                if(c == skipCol) continue;
                for(int r = 0; r < 4; ++r) {
                    if(r == skipRow) continue;
                    m[k++] = cols[static_cast<std::size_t>(c)][r];
                }
            }
            return m[0] * (m[4] * m[8] - m[5] * m[7])
                - m[1] * (m[3] * m[8] - m[5] * m[6])
                + m[2] * (m[3] * m[7] - m[4] * m[6]);
        }
        [[nodiscard]] constexpr T cofactor(int row, int col) const noexcept {
            const T m = minorOf(row, col);
            return ((row + col) & 1) == 0 ? m : -m;
        }
        [[nodiscard]] constexpr T determinant() const noexcept {
            return cols[0][0] * cofactor(0, 0)
                + cols[1][0] * cofactor(0, 1)
                + cols[2][0] * cofactor(0, 2)
                + cols[3][0] * cofactor(0, 3);
        }

        [[nodiscard]] constexpr Mat4 inverted() const noexcept {
            const T det = determinant();
            CAD_VERIFY(det != T(0), "Mat4::inverted on singular matrix");

            Mat4 inv;
            // inverse(col c, row r) = cofactor(row c, col r) / det
            for(int c = 0; c < 4; ++c) {
                for(int r = 0; r < 4; ++r) {
                    inv.cols[static_cast<std::size_t>(c)][r] = cofactor(c, r) / det;
                }
            }
            return inv;
        }

        [[nodiscard]] constexpr bool isIdentity(
            T tol = default_abs_tol<T>) const noexcept
            requires std::is_floating_point_v<T>
        {
            const Mat4<T> i = Identity();
            return almostEqual(cols[0], i.cols[0], tol, tol)
                && almostEqual(cols[1], i.cols[1], tol, tol)
                && almostEqual(cols[2], i.cols[2], tol, tol)
                && almostEqual(cols[3], i.cols[3], tol, tol);
        }

        [[nodiscard]] bool isFinite() const noexcept {
            return cols[0].isFinite() && cols[1].isFinite()
                && cols[2].isFinite() && cols[3].isFinite();
        }
    };

    // aliases
    using Mat4f = Mat4<float>;
    using Mat4d = Mat4<double>;

    // free operators
    template <typename T>
    [[nodiscard]] constexpr Mat4<T> operator*(const Mat4<T>& a, const Mat4<T>& b) noexcept {
        Mat4<T> r;
        for(int i = 0; i < 4; ++i) {
            const Vec4<T>& bc = b.cols[static_cast<std::size_t>(i)];
            r.cols[static_cast<std::size_t>(i)] =
                a.cols[0] * bc.x + a.cols[1] * bc.y + a.cols[2] * bc.z + a.cols[3] * bc.w;
        }
        return r;
    }

    // Full homogeneous transform 
    template <typename T>
    [[nodiscard]] constexpr Vec4<T> operator*(const Mat4<T>& m, const Vec4<T>& v) noexcept {
        return m.cols[0] * v.x + m.cols[1] * v.y + m.cols[2] * v.z + m.cols[3] * v.w;
    }

    template <typename T>
    [[nodiscard]] constexpr Mat4<T> operator*(const Mat4<T>& m, T s) noexcept {
        return Mat4<T>(m.cols[0] * s, m.cols[1] * s, m.cols[2] * s, m.cols[3] * s);
    }

    template <typename T>
    [[nodiscard]] constexpr Mat4<T> operator*(T s, const Mat4<T>& m) noexcept {
        return m * s;
    }

    template <typename T>
    [[nodiscard]] constexpr Mat4<T> operator/(const Mat4<T>& m, T s) noexcept {
        CAD_ASSERT(s != T(0), "Mat4 division by zero");
        return Mat4<T>(m.cols[0] / s, m.cols[1] / s, m.cols[2] / s, m.cols[3] / s);
    }

    // comparison 
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Mat4<T>& a, const Mat4<T>& b) noexcept {
        return a.cols == b.cols;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Mat4<T>& a, const Mat4<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Mat4<T>& a, const Mat4<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.cols[0], b.cols[0], relTol, absTol)
            && almostEqual(a.cols[1], b.cols[1], relTol, absTol)
            && almostEqual(a.cols[2], b.cols[2], relTol, absTol)
            && almostEqual(a.cols[3], b.cols[3], relTol, absTol);
    }

    // free helpers 
    template <typename T>
    [[nodiscard]] constexpr Mat4<T> transpose(const Mat4<T>& m) noexcept {
        return m.transposed();
    }

    template <typename T>
    [[nodiscard]] constexpr Mat4<T> inverse(const Mat4<T>& m) noexcept {
        return m.inverted();
    }

    template <typename T>
    [[nodiscard]] constexpr T determinant(const Mat4<T>& m) noexcept {
        return m.determinant();
    }

    template <typename T>
    [[nodiscard]] constexpr Mat3<T> normalMatrix(const Mat4<T>& m) noexcept {
        return m.upperLeft().inverted().transposed();
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Mat4<T>& m) {
        return os << std::format("{}", m);
    }
} 

namespace std{
    template <typename T>
    struct formatter<cad::Mat4<T>, char> {
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Mat4<T>& m, std::format_context& ctx) const {
            std::string piece;
            piece.reserve(spec.view().size() + 3);
            piece += "{:";
            piece += spec.view();
            piece += '}';
            auto out = ctx.out();
            out = std::format_to(out, "[");
            for(int r = 0; r < 4; ++r) {
                out = std::format_to(out, "{}[", r == 0 ? "" : "\n ");
                for(int c = 0; c < 4; ++c) {
                    if(c != 0) out = std::format_to(out, ", ");
                    out = std::vformat_to(out, piece, std::make_format_args(m.cols[static_cast<std::size_t>(c)][r]));
                }
                out = std::format_to(out, "]");
            }
            return std::format_to(out, "]");
        }
    };
}