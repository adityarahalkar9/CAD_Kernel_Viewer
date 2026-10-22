#pragma once

#include <algorithm>
#include <optional>
#include <type_traits>

#include "../math/Constants.hpp"
#include "../math/Utils.hpp"
#include "../math/Vec2.hpp"
#include "../math/Vec3.hpp"
#include "../math/Vec4.hpp"
#include "../math/Mat3.hpp"
#include "../math/Mat4.hpp"
#include "../math/GeometryTypes.hpp"
#include "../support/Assertions.hpp"

namespace cad {

    // ============================================================================
    //  Result records
    // ============================================================================
    template <typename T>
    struct PointSegmentResult {
        Vec3<T> point{};   // the foot on the segment
        T       t{};       // parameter in [0, 1]: point = a + (b - a) * t
    };

    template <typename T>
    struct PointTriangleResult {
        Vec3<T> point{};   // the foot on the triangle (or its boundary)
        Vec3<T> bary{};    // (w0, w1, w2) w.r.t. (a, b, c); all in [0, 1], sum 1
    };

    template <typename T>
    struct SegmentSegmentResult {
        Vec3<T> point1{};  // foot on segment 1
        Vec3<T> point2{};  // foot on segment 2
        T       s{};       // parameter on segment 1, [0, 1]
        T       t{};       // parameter on segment 2, [0, 1]
        T       distanceSquared{};
    };

    // ============================================================================
    //  Point / segment
    // ============================================================================
    template <typename T>
    [[nodiscard]] PointSegmentResult<T> closestPointOnSegment(const Vec3<T>& p,
                                                              const Vec3<T>& a,
                                                              const Vec3<T>& b) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "closestPointOnSegment requires floating point T");

        const Vec3<T> ab = b - a;
        const T denom = ab.lengthSquared();

        if(denom == T(0)) {
            return {a, T(0)}; // degenerate segment (a == b)
        }

        const T t = std::clamp((p - a).dot(ab) / denom, T(0), T(1));
        return {a + ab * t, t};
    }

    template <typename T>
    [[nodiscard]] T distanceSquaredPointSegment(const Vec3<T>& p, const Vec3<T>& a,
                                                const Vec3<T>& b) noexcept
    {
        const PointSegmentResult<T> r = closestPointOnSegment(p, a, b);
        return (p - r.point).lengthSquared();
    }

    // ============================================================================
    //  Point / triangle - Ericson, region-exact
    //
    //  Classifies p into one of 7 Voronoi regions (3 vertices, 3 edges, interior)
    //  using only dot products - no per-region projections. Exact nearest point
    //  for NON-degenerate triangles; fully degenerate (collinear or collapsed)
    //  inputs still return a sensible foot through the early vertex regions.
    // ============================================================================
    template <typename T>
    [[nodiscard]] PointTriangleResult<T> closestPointOnTriangle(const Vec3<T>& p,
                                                                const Vec3<T>& a,
                                                                const Vec3<T>& b,
                                                                const Vec3<T>& c) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "closestPointOnTriangle requires floating point T");

        const Vec3<T> ab = b - a;
        const Vec3<T> ac = c - a;
        const Vec3<T> ap = p - a;

        const T d1 = ab.dot(ap);
        const T d2 = ac.dot(ap);
        if(d1 <= T(0) && d2 <= T(0)) {
            return {a, Vec3<T>(T(1), T(0), T(0))};              // vertex region A
        }

        const Vec3<T> bp = p - b;
        const T d3 = ab.dot(bp);
        const T d4 = ac.dot(bp);
        if(d3 >= T(0) && d4 <= d3) {
            return {b, Vec3<T>(T(0), T(1), T(0))};              // vertex region B
        }

        const T vc = d1 * d4 - d3 * d2;
        if(vc <= T(0) && d1 >= T(0) && d3 <= T(0)) {
            const T v = d1 / (d1 - d3);
            return {a + ab * v, Vec3<T>(T(1) - v, v, T(0))};    // edge region AB
        }

        const Vec3<T> cp = p - c;
        const T d5 = ab.dot(cp);
        const T d6 = ac.dot(cp);
        if(d6 >= T(0) && d5 <= d6) {
            return {c, Vec3<T>(T(0), T(0), T(1))};              // vertex region C
        }

        const T va = d3 * d6 - d5 * d4;
        if(va <= T(0) && d4 >= T(0) && d6 <= T(0)) {
            const T w = d4 / (d4 - d6);
            return {a + ac * w, Vec3<T>(T(1) - w, T(0), w)};    // edge region AC
        }

        const T vb = d5 * d2 - d1 * d6;
        if(vb <= T(0) && d2 >= T(0) && d6 <= T(0)) {
            const T w = d2 / (d2 - d6);
            return {b + (c - b) * w, Vec3<T>(T(0), T(1) - w, w)}; // edge region BC
        }

        // Interior region.
        const T denom = T(1) / (va + vb + vc);
        const T v = vb * denom;
        const T w = vc * denom;
        return {a + ab * v + ac * w, Vec3<T>(T(1) - v - w, v, w)};
    }

    template <typename T>
    [[nodiscard]] T distanceSquaredPointTriangle(const Vec3<T>& p, const Vec3<T>& a,
                                                 const Vec3<T>& b, const Vec3<T>& c) noexcept
    {
        const PointTriangleResult<T> r = closestPointOnTriangle(p, a, b, c);
        return (p - r.point).lengthSquared();
    }

    // ============================================================================
    //  Point / AABB - per-axis clamp
    // ============================================================================
    template <typename T>
    [[nodiscard]] Vec3<T> closestPointOnAabb(const Vec3<T>& p, const Aabb<T>& box) noexcept
    {
        static_assert(std::is_arithmetic_v<T>,
                      "closestPointOnAabb requires an arithmetic T");
        CAD_ASSERT(!box.isEmpty(), "closestPointOnAabb: box is empty");

        return Vec3<T>(std::clamp(p.x, box.min.x, box.max.x),
                       std::clamp(p.y, box.min.y, box.max.y),
                       std::clamp(p.z, box.min.z, box.max.z));
    }

    template <typename T>
    [[nodiscard]] T distanceSquaredPointAabb(const Vec3<T>& p, const Aabb<T>& box) noexcept
        requires std::is_floating_point_v<T>
    {
        const Vec3<T> q = closestPointOnAabb(p, box);
        return (p - q).lengthSquared();
    }

    // ============================================================================
    //  Point / sphere
    // ============================================================================
    template <typename T>
    [[nodiscard]] Vec3<T> closestPointOnSphere(const Vec3<T>& p,
                                               const Sphere<T>& sphere) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "closestPointOnSphere requires floating point T");

        const Vec3<T> d = p - sphere.center;
        const T len2 = d.lengthSquared();

        if(len2 == T(0)) {
            // Center query: every surface point is equally close; pick +X so the
            // answer is deterministic and documented.
            return sphere.center + Vec3<T>::UnitX() * sphere.radius;
        }

        return sphere.center + d * (sphere.radius / cad::sqrt(len2));
    }

    template <typename T>
    [[nodiscard]] T distancePointSphere(const Vec3<T>& p,
                                        const Sphere<T>& sphere) noexcept
    {
        return std::abs(cad::sqrt((p - sphere.center).lengthSquared()) - sphere.radius);
    }

    // ============================================================================
    //  Point / plane - orthogonal projection
    // ============================================================================
    template <typename T>
    [[nodiscard]] Vec3<T> closestPointOnPlane(const Vec3<T>& p,
                                              const Plane<T>& plane) noexcept
    {
        return plane.projectPoint(p);
    }

    // ============================================================================
    //  Segment / segment - Ericson 5.1.9
    //
    //  Returns the mutually closest point pair, their segment parameters, and
    //  the squared distance. Handles every degeneracy (point segments, parallel
    //  segments) through the clamped-projection fallbacks.
    // ============================================================================
    template <typename T>
    [[nodiscard]] SegmentSegmentResult<T> closestPointsSegmentSegment(
        const Vec3<T>& a1, const Vec3<T>& b1,
        const Vec3<T>& a2, const Vec3<T>& b2) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "closestPointsSegmentSegment requires floating point T");

        const Vec3<T> d1 = b1 - a1;   // direction of segment 1
        const Vec3<T> d2 = b2 - a2;   // direction of segment 2
        const Vec3<T> r = a1 - a2;

        const T a = d1.lengthSquared();
        const T e = d2.lengthSquared();
        const T f = d2.dot(r);

        T s = T(0);
        T t = T(0);

        constexpr T kEps = GeometryTolerance<T>::coincidence;

        if(a <= kEps && e <= kEps) {
            // Both segments are points.
            s = T(0);
            t = T(0);
        }
        else if(a <= kEps) {
            // Segment 1 is a point.
            s = T(0);
            t = std::clamp(f / e, T(0), T(1));
        }
        else {
            const T c = d1.dot(r);
            if(e <= kEps) {
                // Segment 2 is a point.
                t = T(0);
                s = std::clamp(-c / a, T(0), T(1));
            }
            else {
                const T b = d1.dot(d2);
                const T denom = a * e - b * b;

                // Non-parallel segments: closest point of the infinite lines,
                // clamped onto segment 1; parallel (denom ~ 0) picks s = 0.
                s = denom != T(0)
                    ? std::clamp((b * f - c * e) / denom, T(0), T(1))
                    : T(0);

                t = (b * s + f) / e;

                if(t < T(0)) {
                    t = T(0);
                    s = std::clamp(-c / a, T(0), T(1));
                }
                else if(t > T(1)) {
                    t = T(1);
                    s = std::clamp((b - c) / a, T(0), T(1));
                }
            }
        }

        const Vec3<T> p1 = a1 + d1 * s;
        const Vec3<T> p2 = a2 + d2 * t;

        SegmentSegmentResult<T> result;
        result.point1 = p1;
        result.point2 = p2;
        result.s = s;
        result.t = t;
        result.distanceSquared = (p1 - p2).lengthSquared();
        return result;
    }

    template <typename T>
    [[nodiscard]] T distanceSquaredSegmentSegment(const Vec3<T>& a1, const Vec3<T>& b1,
                                                  const Vec3<T>& a2, const Vec3<T>& b2) noexcept
    {
        return closestPointsSegmentSegment(a1, b1, a2, b2).distanceSquared;
    }

} // namespace cad