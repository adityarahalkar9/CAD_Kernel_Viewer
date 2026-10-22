#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <type_traits>
#include <vector>

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
    //  Newell's method
    //
    //  Robust normal estimation for arbitrary simple polygons: convex or concave,
    //  regular or wildly non-convex, exact for planar and least-squares-ish for
    //  near-planar input. Far more stable than cross(firstEdge, secondEdge),
    //  which explodes on collinear first edges.
    //
    //  Returns the UNNORMALIZED area vector (magnitude = 2 * polygon area).
    // ============================================================================
    template <typename T>
    [[nodiscard]] Vec3<T> newellNormal(std::span<const Vec3<T>> polygon) noexcept
    {
        static_assert(std::is_floating_point_v<T>, "newellNormal requires floating point T");
        CAD_VERIFY(polygon.size() >= 3, "newellNormal: polygon needs at least 3 vertices");

        Vec3<T> n = Vec3<T>::Zero();
        const std::size_t count = polygon.size();

        for(std::size_t i = 0; i < count; ++i) {
            const Vec3<T>& v = polygon[i];
            const Vec3<T>& vnext = polygon[(i + 1) % count];

            n.x += (v.y - vnext.y) * (v.z + vnext.z);
            n.y += (v.z - vnext.z) * (v.x + vnext.x);
            n.z += (v.x - vnext.x) * (v.y + vnext.y);
        }

        return n;
    }

    template <typename T>
    [[nodiscard]] Vec3<T> polygonNormal(std::span<const Vec3<T>> polygon) noexcept
    {
        const Vec3<T> n = newellNormal(polygon);
        CAD_ASSERT(n.lengthSquared() > T(0),
                   "polygonNormal: polygon has zero area (degenerate or collinear)");
        return n.normalized();
    }

    template <typename T>
    [[nodiscard]] T polygonArea(std::span<const Vec3<T>> polygon) noexcept
    {
        return newellNormal(polygon).length() / T(2);
    }

    // Best-fit plane: Newell normal + plane through the vertex average.
    // Exact for planar polygons; for near-planar input this is the natural
    // least-squares-ish compromise (documented, deterministic).
    template <typename T>
    [[nodiscard]] Plane<T> polygonPlane(std::span<const Vec3<T>> polygon) noexcept
    {
        const Vec3<T> n = polygonNormal(polygon);

        Vec3<T> mean = Vec3<T>::Zero();
        for(const Vec3<T>& v : polygon) {
            mean += v;
        }
        mean = mean / static_cast<T>(polygon.size());

        return Plane<T>::fromPointAndNormal(mean, n);
    }

    // Area-weighted centroid via fan triangulation with SIGNED areas measured
    // against the polygon's own Newell normal. Exact for any SIMPLE (concave
    // included) planar polygon: overlapping fan pieces cancel through the signs.
    template <typename T>
    [[nodiscard]] Vec3<T> polygonCentroid(std::span<const Vec3<T>> polygon) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "polygonCentroid requires floating point T");
        CAD_VERIFY(polygon.size() >= 3, "polygonCentroid: needs at least 3 vertices");

        const Vec3<T> nHat = polygonNormal(polygon);

        T   total2A = T(0);          // accumulates 2 * total signed area
        Vec3<T> acc = Vec3<T>::Zero();

        const Vec3<T>& v0 = polygon[0];
        for(std::size_t i = 1; i + 1 < polygon.size(); ++i) {
            const Vec3<T>& a = v0;
            const Vec3<T>& b = polygon[i];
            const Vec3<T>& c = polygon[i + 1];

            const T signed2A = (b - a).cross(c - a).dot(nHat); // 2 * signed area

            total2A += signed2A;
            acc += (a + b + c) * signed2A;   // 3 * centroid * (2 * signed area)
        }

        CAD_VERIFY(total2A != T(0),
                   "polygonCentroid: polygon has zero signed area (degenerate)");

        return acc / (T(3) * total2A);
    }

    // ============================================================================
    //  Fan triangulation
    //
    //  Emits index triples (0, i, i+1) over the polygon's vertex order. Only
    //  produces a valid triangulation for CONVEX polygons; concave faces need a
    //  real ear-clipper later (roadmap item). Still useful everywhere convex
    //  data is guaranteed: box/sphere/cylinder caps, clipped cross-sections.
    // ============================================================================
    template <typename Index>
    [[nodiscard]] std::vector<std::array<Index, 3>> triangulateFan(
        std::size_t vertexCount)
    {
        static_assert(std::is_integral_v<Index>, "triangulateFan requires an integer Index");
        CAD_VERIFY(vertexCount >= 3, "triangulateFan: needs at least 3 vertices");

        std::vector<std::array<Index, 3>> triangles;
        triangles.reserve(vertexCount - 2);
        for(std::size_t i = 1; i + 1 < vertexCount; ++i) {
            triangles.push_back({Index(0), Index(i), Index(i + 1)});
        }
        return triangles;
    }

    // ============================================================================
    //  Point / plane classification
    // ============================================================================
    enum class PlaneSide : std::uint8_t {
        Negative = 0,   // opposite side from the normal
        On = 1,
        Positive = 2,   // the side the normal points to
    };

    template <typename T>
    [[nodiscard]] PlaneSide classifyPointToPlane(
        const Plane<T>& plane, const Vec3<T>& p,
        T tolerance = GeometryTolerance<T>::point_on_plane) noexcept
    {
        const T sd = plane.signedDistance(p);
        if(sd > tolerance) return PlaneSide::Positive;
        if(sd < -tolerance) return PlaneSide::Negative;
        return PlaneSide::On;
    }

    // ============================================================================
    //  Segment / plane intersection