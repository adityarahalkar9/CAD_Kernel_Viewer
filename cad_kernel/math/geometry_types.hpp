#pragma once

#include <array>
#include <cstddef>
#include <format>
#include <ostream>
#include <span>
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
#include "../support/Assertions.hpp"

namespace cad{
    // Ray
    template <typename T>
    struct Ray {
        static_assert(std::is_floating_point_v<T>, "Ray<T> requires a floating point T");

        Vec3<T> origin{};
        Vec3<T> direction{};       // must be unit length
        T       tmin{T(0)};        // inclusive lower hit bound (self-hit epsilon)
        T       tmax{infinity<T>()}; // inclusive upper hit bound

        constexpr Ray() noexcept = default;

        constexpr Ray(const Vec3<T>& origin, const Vec3<T>& unitDirection) noexcept
            : origin(origin), direction(unitDirection) {
        }

        constexpr Ray(const Vec3<T>& origin, const Vec3<T>& unitDirection,
                      T tmin, T tmax) noexcept
            : origin(origin), direction(unitDirection), tmin(tmin), tmax(tmax) {
        }

        [[nodiscard]] constexpr Vec3<T> pointAt(T t) const noexcept {
            return origin + direction * t;
        }

        // Loose validity: finite, direction ~unit, ordered interval.
        [[nodiscard]] bool isValid() const noexcept {
            return origin.isFinite()
                && direction.isFinite()
                && almostZero(direction.lengthSquared() - T(1), T(1e-4))
                && tmax >= tmin;
        }
    };

    //  Plane:  dot(normal, p) == distance   for every point p on the plane
    template <typename T>
    struct Plane {
        static_assert(std::is_floating_point_v<T>, "Plane<T> requires a floating point T");

        Vec3<T> normal{};   // expected unit length
        T       distance{T(0)};

        constexpr Plane() noexcept = default;

        constexpr Plane(const Vec3<T>& unitNormal, T distance) noexcept
            : normal(unitNormal), distance(distance) {
        }

        [[nodiscard]] static Plane fromPointAndNormal(const Vec3<T>& point,
                                                      const Vec3<T>& unitNormal) noexcept {
            return Plane(unitNormal, unitNormal.dot(point));
        }

        // Plane through three points. CCW winding (a, b, c) viewed from the
        // normal side - the same convention as mesh triangle face.
        // Degenerate (collinear) input trips CAD_VERIFY.
        [[nodiscard]] static Plane fromThreePoints(const Vec3<T>& a, const Vec3<T>& b,
                                                   const Vec3<T>& c) noexcept {
            const Vec3<T> n = b.cross(c); // placeholder, replaced below
            (void)n;
            const Vec3<T> u = b - a;
            const Vec3<T> v = c - a;
            const Vec3<T> cross = u.cross(v);
            const T       len2 = cross.lengthSquared();
            CAD_VERIFY(len2 > T(0), "Plane::fromThreePoints: collinear points");
            const Vec3<T> unit = cross / cad::sqrt(len2);
            return Plane(unit, unit.dot(a));
        }

        // Positive halfspace: the side the normal points to.
        [[nodiscard]] constexpr T signedDistance(const Vec3<T>& p) const noexcept {
            return normal.dot(p) - distance;
        }

        // Foot of p on the plane.
        [[nodiscard]] constexpr Vec3<T> projectPoint(const Vec3<T>& p) const noexcept {
            return p - normal * signedDistance(p);
        }

        // Closest point to the coordinate origin.
        [[nodiscard]] constexpr Vec3<T> origin() const noexcept {
            return normal * distance;
        }

        // Same plane, opposite orientation.
        [[nodiscard]] constexpr Plane flipped() const noexcept {
            return Plane(-normal, -distance);
        }

        // Normalize a non-unit plane (scales distance by 1/|n|).
        [[nodiscard]] Plane normalized() const noexcept {
            const T len = normal.length();
            CAD_ASSERT(len > T(0), "Plane::normalized on zero normal");
            return Plane(normal / len, distance / len);
        }

        [[nodiscard]] constexpr bool isNormalized(
            T tol = GeometryTolerance<T>::unit_deviation) const noexcept {
            return normal.isNormalized(tol);
        }

        [[nodiscard]] bool isValid() const noexcept {
            return normal.isFinite() && std::isfinite(distance)
                && normal.lengthSquared() > T(0);
        }
    };

    // Aabb - axis-aligned bounding box
    template <typename T>
    struct Aabb {
        static_assert(std::is_arithmetic_v<T>, "Aabb<T> requires an arithmetic T");

        Vec3<T> min{};
        Vec3<T> max{};

        constexpr Aabb() noexcept = default;

        constexpr Aabb(const Vec3<T>& minCorner, const Vec3<T>& maxCorner) noexcept
            : min(minCorner), max(maxCorner) {
        }

        // The identity of union: a box containing no points.
        [[nodiscard]] static Aabb Empty() noexcept
            requires std::is_floating_point_v<T>
        {
            return Aabb(Vec3<T>(infinity<T>(), infinity<T>(), infinity<T>()),
                        Vec3<T>(-infinity<T>(), -infinity<T>(), -infinity<T>()));
        }

        [[nodiscard]] static Aabb FromPoint(const Vec3<T>& p) noexcept {
            return Aabb(p, p);
        }

        // Two arbitrary corners; order is fixed internally.
        [[nodiscard]] static Aabb FromPointPair(const Vec3<T>& a, const Vec3<T>& b) noexcept {
            return Aabb(a.cwMin(b), a.cwMax(b));
        }

        [[nodiscard]] static Aabb FromCenterHalfExtents(const Vec3<T>& center,
                                                        const Vec3<T>& halfExtents) noexcept {
            return Aabb(center - halfExtents, center + halfExtents);
        }

        // Envelope of a point set.
        [[nodiscard]] static Aabb fromPoints(std::span<const Vec3<T>> points) noexcept
            requires std::is_floating_point_v<T>
        {
            Aabb box = Empty();
            for(const Vec3<T>& p : points) {
                box.extend(p);
            }
            return box;
        }

        // -- growing (in place) ---------------------------------------------------------
        constexpr Aabb& extend(const Vec3<T>& p) noexcept {
            min = min.cwMin(p);
            max = max.cwMax(p);
            return *this;
        }

        constexpr Aabb& expand(const Aabb& other) noexcept {
            min = min.cwMin(other.min);
            max = max.cwMax(other.max);
            return *this;
        }

        [[nodiscard]] constexpr Aabb extended(const Vec3<T>& p) const noexcept {
            return Aabb(min.cwMin(p), max.cwMax(p));
        }

        [[nodiscard]] constexpr Aabb mergedWith(const Aabb& other) const noexcept {
            return Aabb(min.cwMin(other.min), max.cwMax(other.max));
        }

        // -- queries ----------------------------------------------------------------------
        [[nodiscard]] constexpr bool isEmpty() const noexcept {
            return min.x > max.x || min.y > max.y || min.z > max.z;
        }

        // Finite corners in the canonical order min <= max (Empty() is invalid by
        // construction - callers building boxes through extend() never see that).
        [[nodiscard]] constexpr bool isValid() const noexcept
            requires std::is_floating_point_v<T>
        {
            return min.isFinite() && max.isFinite()
                && min.x <= max.x && min.y <= max.y && min.z <= max.z;
        }

        [[nodiscard]] constexpr Vec3<T> center() const noexcept
            requires std::is_floating_point_v<T>
        {
            return (min + max) * T(0.5);
        }

        // full size per axis
        [[nodiscard]] constexpr Vec3<T> extents() const noexcept {
            return max - min;
        }

        [[nodiscard]] constexpr Vec3<T> halfExtents() const noexcept
            requires std::is_floating_point_v<T>
        {
            return extents() * T(0.5);
        }

        [[nodiscard]] constexpr Vec3<T> diagonal() const noexcept {
            return max - min;
        }

        [[nodiscard]] constexpr T volume() const noexcept
            requires std::is_floating_point_v<T>
        {
            const Vec3<T> e = extents();
            return (e.x < T(0) || e.y < T(0) || e.z < T(0)) ? T(0) : e.x * e.y * e.z;
        }

        // Surface area (SAH cost primitive for the AABB tree in Part 8).
        [[nodiscard]] constexpr T surfaceArea() const noexcept
            requires std::is_floating_point_v<T>
        {
            const Vec3<T> e = extents();
            if(e.x < T(0) || e.y < T(0) || e.z < T(0)) {
                return T(0);
            }
            return T(2) * (e.x * e.y + e.y * e.z + e.z * e.x);
        }

        // Axis with the largest extent (split-axis heuristic); ties -> lowest axis.
        [[nodiscard]] constexpr int longestAxis() const noexcept {
            const Vec3<T> e = extents();
            if(e.x >= e.y && e.x >= e.z) return 0;
            if(e.y >= e.z)               return 1;
            return 2;
        }

        [[nodiscard]] constexpr T maxExtent() const noexcept {
            return extents().maxComponent();
        }

        // Inclusive containment.
        [[nodiscard]] constexpr bool containsPoint(const Vec3<T>& p) const noexcept {
            return p.x >= min.x && p.x <= max.x
                && p.y >= min.y && p.y <= max.y
                && p.z >= min.z && p.z <= max.z;
        }

        // Inclusive: touching boxes count as intersecting (closed-set semantics,
        // which keeps BVH traversal sound).
        [[nodiscard]] constexpr bool intersectsBox(const Aabb& other) const noexcept {
            return min.x <= other.max.x && max.x >= other.min.x
                && min.y <= other.max.y && max.y >= other.min.y
                && min.z <= other.max.z && max.z >= other.min.z;
        }

        [[nodiscard]] constexpr bool containsBox(const Aabb& other) const noexcept {
            return min.x <= other.min.x && max.x >= other.max.x
                && min.y <= other.min.y && max.y >= other.max.y
                && min.z <= other.min.z && max.z >= other.max.z;
        }

        // derived boxes
        [[nodiscard]] constexpr Aabb translated(const Vec3<T>& offset) const noexcept {
            return Aabb(min + offset, max + offset);
        }

        // Uniform inflation by an epsilon (BVH robustness padding).
        [[nodiscard]] constexpr Aabb padded(T epsilon) const noexcept {
            return Aabb(min - Vec3<T>(epsilon, epsilon, epsilon),
                        max + Vec3<T>(epsilon, epsilon, epsilon));
        }

        // Conservative bound of this box under an affine transform: transforms the
        // 8 corners, takes their envelope. EXACT for translations/axis-aligned
        // scales; conservative (>= exact rotated box) otherwise.
        [[nodiscard]] Aabb transformed(const Mat4<T>& m) const noexcept
            requires std::is_floating_point_v<T>
        {
            Aabb box = Empty();
            for(const Vec3<T>& c : corners()) {
                box.extend(m.transformPoint(c));
            }
            return box;
        }

        // Corner order: bit i selects max over min on axis i
        //   index & 1 -> x, index & 2 -> y, index & 4 -> z
        [[nodiscard]] constexpr std::array<Vec3<T>, 8> corners() const noexcept {
            return {
                Vec3<T>(min.x, min.y, min.z),
                Vec3<T>(max.x, min.y, min.z),
                Vec3<T>(min.x, max.y, min.z),
                Vec3<T>(max.x, max.y, min.z),
                Vec3<T>(min.x, min.y, max.z),
                Vec3<T>(max.x, min.y, max.z),
                Vec3<T>(min.x, max.y, max.z),
                Vec3<T>(max.x, max.y, max.z),
            };
        }
    };

    // Sphere
    template <typename T>
    struct Sphere {
        static_assert(std::is_floating_point_v<T>, "Sphere<T> requires a floating point T");

        Vec3<T> center{};
        T       radius{T(0)};

        constexpr Sphere() noexcept = default;

        constexpr Sphere(const Vec3<T>& center, T radius) noexcept
            : center(center), radius(radius) {
        }

        [[nodiscard]] constexpr bool containsPoint(const Vec3<T>& p) const noexcept {
            return distanceSquared(center, p) <= radius * radius;
        }

        [[nodiscard]] constexpr T volume() const noexcept {
            return (T(4) / T(3)) * pi_v<T> *radius * radius * radius;
        }

        [[nodiscard]] constexpr T surfaceArea() const noexcept {
            return T(4) * pi_v<T> *radius * radius;
        }

        // Exact minimal sphere enclosing both inputs (two-sphere problem).
        // Touching/contained cases keep the larger sphere untouched.
        [[nodiscard]] static Sphere merged(const Sphere& a, const Sphere& b) noexcept {
            const Vec3<T> d = b.center - a.center;
            const T       dist2 = d.lengthSquared();

            if(dist2 <= (a.radius - b.radius) * (a.radius - b.radius)) {
                return a; // b inside a
            }
            if(dist2 <= (b.radius - a.radius) * (b.radius - a.radius)) {
                return b; // a inside b
            }

            const T dist = cad::sqrt(dist2);
            const T newR = (dist + a.radius + b.radius) / T(2);
            const Vec3<T> dir = dist > T(0) ? d / dist : Vec3<T>::UnitX();
            return Sphere(a.center + dir * (newR - a.radius), newR);
        }

        // Loose growth to include a point (center kept) - quick bound updates.
        [[nodiscard]] constexpr Sphere extended(const Vec3<T>& p) const noexcept {
            const T d = cad::sqrt(distanceSquared(center, p));
            return d <= radius ? *this : Sphere(center, d);
        }
    };

    // aliases
    using Rayd = Ray<double>;   using Rayf = Ray<float>;
    using Planed = Plane<double>; using Planef = Plane<float>;
    using Aabbd = Aabb<double>;  using Aabbf = Aabb<float>;
    using Aabbi = Aabb<std::int32_t>;
    using Sphered = Sphere<double>; using Spheref = Sphere<float>;

    // comparisons 
    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Ray<T>& a, const Ray<T>& b) noexcept {
        return a.origin == b.origin && a.direction == b.direction
            && a.tmin == b.tmin && a.tmax == b.tmax;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Ray<T>& a, const Ray<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Plane<T>& a, const Plane<T>& b) noexcept {
        return a.normal == b.normal && a.distance == b.distance;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Plane<T>& a, const Plane<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Aabb<T>& a, const Aabb<T>& b) noexcept {
        return a.min == b.min && a.max == b.max;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Aabb<T>& a, const Aabb<T>& b) noexcept {
        return !(a == b);
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator==(const Sphere<T>& a, const Sphere<T>& b) noexcept {
        return a.center == b.center && a.radius == b.radius;
    }

    template <typename T>
    [[nodiscard]] constexpr bool operator!=(const Sphere<T>& a, const Sphere<T>& b) noexcept {
        return !(a == b);
    }

    // approximate equality 
    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Ray<T>& a, const Ray<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        return almostEqual(a.origin, b.origin, relTol, absTol)
            && almostEqual(a.direction, b.direction, relTol, absTol)
            && almostEqual(a.tmin, b.tmin, relTol, absTol)
            && almostEqual(a.tmax, b.tmax, relTol, absTol);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Plane<T>& a, const Plane<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        return almostEqual(a.normal, b.normal, relTol, absTol)
            && almostEqual(a.distance, b.distance, relTol, absTol);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Aabb<T>& a, const Aabb<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept
        requires std::is_floating_point_v<T>
    {
        return almostEqual(a.min, b.min, relTol, absTol)
            && almostEqual(a.max, b.max, relTol, absTol);
    }

    template <typename T>
    [[nodiscard]] constexpr bool almostEqual(const Sphere<T>& a, const Sphere<T>& b,
                                             T relTol = default_rel_tol<T>,
                                             T absTol = default_abs_tol<T>) noexcept {
        return almostEqual(a.center, b.center, relTol, absTol)
            && almostEqual(a.radius, b.radius, relTol, absTol);
    }

    // -- stream output -----------------------------------------------------------------------
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Ray<T>& r) {
        return os << std::format("{}", r);
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Plane<T>& p) {
        return os << std::format("{}", p);
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Aabb<T>& b) {
        return os << std::format("{}", b);
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& os, const Sphere<T>& s) {
        return os << std::format("{}", s);
    }

}

namespace std {

    template <typename T>
    struct formatter<cad::Ray<T>, char> {
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Ray<T>& r, std::format_context& ctx) const {
            std::string piece;
            piece += "{:"; piece += spec.view(); piece += '}';
            const std::string text = std::vformat(
                "[o " + piece + ", d " + piece + ", t [" + piece + ", " + piece + "]]",
                std::make_format_args(r.origin, r.direction, r.tmin, r.tmax));
            return std::format_to(ctx.out(), "{}", text);
        }
    };

    template <typename T>
    struct formatter<cad::Plane<T>, char>{
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx){
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Plane<T>& p, std::format_context& ctx) const{
            std::string piece;
            piece += "{:"; piece += spec.view(); piece += '}';
            const std::string text = std::vformat(
                "[n " + piece + ", d " + piece + "]",
                std::make_format_args(p.normal, p.distance));
            return std::format_to(ctx.out(), "{}", text);
        }
    };

    template <typename T>
    struct formatter<cad::Aabb<T>, char>{
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Aabb<T>& b, std::format_context& ctx) const {
            std::string piece;
            piece += "{:"; piece += spec.view(); piece += '}';
            const std::string text = std::vformat(
                "[min " + piece + ", max " + piece + "]",
                std::make_format_args(b.min, b.max));
            return std::format_to(ctx.out(), "{}", text);
        }
    };

    template <typename T>
    struct formatter<cad::Sphere<T>, char> {
        cad::detail::FormatSpec spec;
        constexpr auto parse(std::format_parse_context& ctx) {
            return cad::detail::parseSpec(spec, ctx);
        }
        auto format(const cad::Sphere<T>& s, std::format_context& ctx) const {
            std::string piece;
            piece += "{:"; piece += spec.view(); piece += '}';
            const std::string text = std::vformat(
                "[c " + piece + ", r " + piece + "]",
                std::make_format_args(s.center, s.radius));
            return std::format_to(ctx.out(), "{}", text);
        }
    };

} // namespace std