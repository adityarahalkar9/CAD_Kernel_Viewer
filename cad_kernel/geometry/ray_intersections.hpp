#pragma once

#include <algorithm>
#include <optional>
#include <type_traits>
#include <utility>

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
    //  Hit record
    // ============================================================================
    template <typename T>
    struct RayHit {
        static_assert(std::is_floating_point_v<T>, "RayHit<T> requires floating point T");

        T       t{};          // ray parameter of the hit
        Vec3<T> point{};      // origin + direction * t
        Vec3<T> normal{};     // unit geometric normal (triangle: winding order;
        // sphere: outward; plane: plane normal) - NOT flipped
        Vec3<T> bary{};       // triangle only: (w0, w1, w2) for (v0, v1, v2); else zeros
        bool    backface{false}; // true when the ray arrives from behind the face
        // (sphere: hit from inside; plane: denom > 0)
    };

    using RayHitd = RayHit<double>;
    using RayHitf = RayHit<float>;

    // ============================================================================
    //  Watertight ray frame
    //
    //  Precomputes the shear that projects the ray onto the plane perpendicular
    //  to its DOMINANT direction axis. In sheared space the ray is the +W axis
    //  through the origin, and containment reduces to 2D edge functions.
    //
    //  Build ONCE per ray and test every triangle against it: that shared frame
    //  (identical rounding for identical vertices) is half the watertightness
    //  guarantee; the antisymmetric edge evaluation below is the other half.
    // ============================================================================
    template <typename T>
    struct WatertightRay {
        static_assert(std::is_floating_point_v<T>,
                      "WatertightRay<T> requires floating point T");

        Vec3<T> origin{};     // ray origin (unsheared; kept for hit assembly)
        Vec3<T> direction{};  // unit direction (kept for hit assembly + backface)

        int k{2};             // dominant axis of |direction|
        int iu{0};            // (k + 1) % 3 - sheared U axis
        int iv{1};            // (k + 2) % 3 - sheared V axis

        T Sx{};               // direction[iu] / direction[k]
        T Sy{};               // direction[iv] / direction[k]
        T Sz{};               // 1 / direction[k]

        T tmin{T(0)};         // active interval (copied from the ray)
        T tmax{infinity<T>()};

        constexpr WatertightRay() noexcept = default;

        explicit WatertightRay(const Ray<T>& ray) noexcept
            : origin(ray.origin), direction(ray.direction),
            tmin(ray.tmin), tmax(ray.tmax)
        {
            CAD_ASSERT(ray.direction.isNormalized(T(1e-3)),
                       "WatertightRay: direction must be unit length");
            CAD_ASSERT(ray.direction.isFinite(), "WatertightRay: non-finite direction");

            const Vec3<T> a = ray.direction.cwAbs();
            k = a.maxAxis();          // |d[k]| >= 1/sqrt(3) for any unit vector
            iu = (k + 1) % 3;
            iv = (k + 2) % 3;

            CAD_ASSERT(ray.direction[k] != T(0),
                       "WatertightRay: dominant direction component is zero");

            Sz = T(1) / ray.direction[k];
            Sx = ray.direction[iu] * Sz;
            Sy = ray.direction[iv] * Sz;
        }
    };

    namespace detail {

        // One vertex, sheared into the ray frame. U/V are the projected 2D
        // coordinates, W is the depth along the (unnormalized) ray axis.
        template <typename T>
        struct ShearedVertex {
            T u{};
            T v{};
            T w{};
        };

        template <typename T>
        [[nodiscard]] inline ShearedVertex<T> shearVertex(const WatertightRay<T>& wr,
                                                          const Vec3<T>& vertex) noexcept
        {
            const Vec3<T> rel = vertex - wr.origin;
            const T w = rel[wr.k];
            return ShearedVertex<T>{
                rel[wr.iu] - wr.Sx * w,
                    rel[wr.iv] - wr.Sy * w,
                    w};
        }

    } // namespace detail

    // ============================================================================
    //  Watertight ray/triangle intersection
    //
    //  Returns nullopt for misses (including degenerate triangles: det == 0
    //  swallows collinear vertices AND rays parallel to the triangle plane
    //  without ever producing NaN).
    //
    //  cullBackfaces: reject hits arriving from behind the winding normal
    //  (front = dot(direction, geometricNormal) < 0).
    //
    //  THE EDGE FUNCTION LINE IS LOAD-BEARING:
    //
    //      const T e = a.u * b.v - a.v * b.u;
    //
    //  * NO std::fma here, ever. Each product must round SEPARATELY so that
    //    cross2(P, Q) and cross2(Q, P) are exact negations of each other for
    //    bit-identical inputs (IEEE multiplication is commutative and
    //    negation-symmetric; subtraction is antisymmetric). An FMA-fused
    //    version breaks that pairing and reopens cracks along shared edges.
    //  * The sheared coordinates of a shared vertex are bit-identical in every
    //    incident triangle (same WatertightRay, same operations), so the
    //    negation pairing actually fires across triangle boundaries.
    // ============================================================================
    template <typename T>
    [[nodiscard]] std::optional<RayHit<T>> intersectTriangleWatertight(
        const WatertightRay<T>& wr,
        const Vec3<T>& v0, const Vec3<T>& v1, const Vec3<T>& v2,
        bool cullBackfaces = false) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "intersectTriangleWatertight requires floating point T");

        // -- shear all vertices into the shared ray frame ---------------------------
        const detail::ShearedVertex<T> a = detail::shearVertex<T>(wr, v0);
        const detail::ShearedVertex<T> b = detail::shearVertex<T>(wr, v1);
        const detail::ShearedVertex<T> c = detail::shearVertex<T>(wr, v2);

        // -- edge functions (2D crosses; see the load-bearing note above) -----------
        const T e_ab = a.u * b.v - a.v * b.u;   // edge A->B (opposite C)
        const T e_bc = b.u * c.v - b.v * c.u;   // edge B->C (opposite A)
        const T e_ca = c.u * a.v - c.v * a.u;   // edge C->A (opposite B)

        const T det = e_ab + e_bc + e_ca;

        // Degenerate triangle or ray parallel to the triangle plane.
        // NaN inputs would also land here (NaN != 0 is true, so they do NOT -
        // callers CAD_ASSERT finite inputs upstream; the BVH does).
        if(det == T(0)) {
            return std::nullopt;
        }

        // Containment: all edge values must share a sign (zeros allowed - they
        // are exact edge/vertex hits, the watertight cases we care about).
        const T minE = std::min({e_ab, e_bc, e_ca});
        const T maxE = std::max({e_ab, e_bc, e_ca});
        if(minE < T(0) && maxE > T(0)) {
            return std::nullopt;
        }

        // -- backface classification ------------------------------------------------
        // sign(det) == sign(dot(direction, normal)) * sign(direction[k]).
        const bool backface = (det > T(0)) == (wr.direction[wr.k] > T(0));
        if(cullBackfaces && backface) {
            return std::nullopt;
        }

        // -- hit depth --------------------------------------------------------------
        // t = (interpolated sheared depth) / direction[k]; barycentric weights are
        // the edge values normalized by det (edge opposite a vertex weights it).
        const T t = ((a.w * e_bc + b.w * e_ca + c.w * e_ab) / det) * wr.Sz;

        if(t < wr.tmin || t > wr.tmax) {
            return std::nullopt;
        }

        // -- assemble the hit ---------------------------------------------------------
        const Vec3<T> n = (v1 - v0).cross(v2 - v0).normalized();
        // det != 0 implies the vertices are not collinear, hence n is well
        // defined (the sheared projection has nonzero area <=> non-collinear).

        RayHit<T> hit;
        hit.t = t;
        hit.point = wr.origin + wr.direction * t;
        hit.normal = n;
        hit.bary = Vec3<T>(e_bc / det, e_ca / det, e_ab / det);
        hit.backface = backface;
        return hit;
    }

    // Convenience overload for one-off queries (tests, tools). Hot paths (the
    // AABB tree in Part 8) reuse a precomputed WatertightRay across the whole
    // traversal instead.
    template <typename T>
    [[nodiscard]] std::optional<RayHit<T>> intersectTriangleWatertight(
        const Ray<T>& ray,
        const Vec3<T>& v0, const Vec3<T>& v1, const Vec3<T>& v2,
        bool cullBackfaces = false) noexcept
    {
        const WatertightRay<T> wr(ray);
        return intersectTriangleWatertight(wr, v0, v1, v2, cullBackfaces);
    }

    // ============================================================================
    //  Moller-Trumbore intersection (classic)
    //
    //  Kept for cross-validation: Part 14's test suite throws the same queries at
    //  both algorithms and compares. It is FINE for interior hits - it is only
    //  unreliable in the exact-shared-edge regime, which is exactly where the
    //  watertight version earns its keep.
    // ============================================================================
    template <typename T>
    [[nodiscard]] std::optional<RayHit<T>> intersectTriangleMollerTrumbore(
        const Ray<T>& ray,
        const Vec3<T>& v0, const Vec3<T>& v1, const Vec3<T>& v2,
        T epsilon = GeometryTolerance<T>::parallel) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "intersectTriangleMollerTrumbore requires floating point T");

        const Vec3<T> e1 = v1 - v0;
        const Vec3<T> e2 = v2 - v0;

        const Vec3<T> pvec = ray.direction.cross(e2);
        const T det = e1.dot(pvec);

        if(det > -epsilon && det < epsilon) {
            return std::nullopt; // parallel, or degenerate triangle
        }

        const T invDet = T(1) / det;

        const Vec3<T> tvec = ray.origin - v0;
        const T u = tvec.dot(pvec) * invDet;
        if(u < T(0) || u > T(1)) {
            return std::nullopt;
        }

        const Vec3<T> qvec = tvec.cross(e1);
        const T v = ray.direction.dot(qvec) * invDet;
        if(v < T(0) || u + v > T(1)) {
            return std::nullopt;
        }

        const T t = e2.dot(qvec) * invDet;
        if(t < ray.tmin || t > ray.tmax) {
            return std::nullopt;
        }

        RayHit<T> hit;
        hit.t = t;
        hit.point = ray.origin + ray.direction * t;
        hit.normal = e1.cross(e2).normalized();
        hit.bary = Vec3<T>(T(1) - u - v, u, v);
        hit.backface = det > T(0); // det == dot(direction, geometricNormal)
        return hit;
    }

    // ============================================================================
    //  Ray / AABB - slab test, NaN-proof variant
    //
    //  The textbook "branchless" slab test divides by the inverse direction and
    //  produces NaN whenever a slab boundary difference is 0 and the component
    //  is 0 (0 * inf). std::min/std::max do NOT discard NaN, so a NaN silently
    //  poisons the interval and boxes get missed. This version only divides on
    //  axes with a nonzero component and handles parallel axes explicitly - a
    //  handful of well-predicted branches, in exchange for never lying.
    //
    //  Contract:
    //   * returns true when [t0Out, t1Out] (the slab interval, clamped to the
    //     ray's active [tmin, tmax]) is non-empty.
    //   * an origin INSIDE the box reports a hit (interval starts at ray.tmin).
    //   * segment queries: pass a ray with tmin = 0, tmax = segment length.
    // ============================================================================
    template <typename T>
    [[nodiscard]] bool intersectAabb(const Ray<T>& ray, const Aabb<T>& box,
                                     T& t0Out, T& t1Out) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "intersectAabb requires floating point T");

        T tNear = -infinity<T>();
        T tFar = infinity<T>();

        for(int axis = 0; axis < 3; ++axis) {
            const T o = ray.origin[axis];
            const T d = ray.direction[axis];
            const T lo = box.min[axis];
            const T hi = box.max[axis];

            if(d != T(0)) {
                // Division only on non-degenerate axes: no 0/0, no 0*inf.
                // Overflowing quotients become +/-inf, which the interval
                // arithmetic below handles correctly.
                T s0 = (lo - o) / d;
                T s1 = (hi - o) / d;
                if(s0 > s1) {
                    std::swap(s0, s1);
                }
                tNear = std::max(tNear, s0);
                tFar = std::min(tFar, s1);
            }
            else {
                // Parallel to this slab: hit only if the origin is inside it.
                if(o < lo || o > hi) {
                    return false;
                }
            }
        }

        t0Out = std::max(tNear, ray.tmin);
        t1Out = std::min(tFar, ray.tmax);
        return t0Out <= t1Out;
    }

    template <typename T>
    [[nodiscard]] bool intersectAabb(const Ray<T>& ray, const Aabb<T>& box) noexcept
    {
        T t0, t1;
        return intersectAabb(ray, box, t0, t1);
    }

    // ============================================================================
    //  Ray / sphere - geometric quadratic form
    //
    //  Solves |origin + t*direction - center|^2 = r^2 with direction unit:
    //  t^2 + 2 b t + c = 0 where b = dot(origin - center, direction).
    //  The reduced form avoids computing `a` (== 1) entirely.
    // ============================================================================
    template <typename T>
    [[nodiscard]] std::optional<RayHit<T>> intersectSphere(const Ray<T>& ray,
                                                           const Sphere<T>& sphere) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "intersectSphere requires floating point T");
        CAD_ASSERT(sphere.radius >= T(0), "intersectSphere: negative radius");

        const Vec3<T> oc = ray.origin - sphere.center;
        const T b = oc.dot(ray.direction);
        const T c = oc.lengthSquared() - sphere.radius * sphere.radius;

        const T disc = b * b - c;
        if(disc < T(0)) {
            return std::nullopt;
        }

        const T s = cad::sqrt(disc);
        const T t0 = -b - s;   // near root
        const T t1 = -b + s;   // far root

        // Prefer the near root inside the active interval; an origin inside the
        // sphere naturally falls through to the far (exit) root.
        const T t = (t0 >= ray.tmin) ? t0 : t1;
        if(t < ray.tmin || t > ray.tmax) {
            return std::nullopt;
        }

        RayHit<T> hit;
        hit.t = t;
        hit.point = ray.origin + ray.direction * t;
        hit.normal = (hit.point - sphere.center) / sphere.radius;
        hit.backface = c < T(0); // started inside the sphere
        return hit;
    }

    // ============================================================================
    //  Ray / plane
    // ============================================================================
    template <typename T>
    [[nodiscard]] std::optional<RayHit<T>> intersectPlane(const Ray<T>& ray,
                                                          const Plane<T>& plane) noexcept
    {
        static_assert(std::is_floating_point_v<T>,
                      "intersectPlane requires floating point T");

        const T denom = plane.normal.dot(ray.direction);

        // Parallel (or anti-parallel): no unique intersection point.
        if(denom > -GeometryTolerance<T>::parallel &&
           denom < GeometryTolerance<T>::parallel) {
            return std::nullopt;
        }

        const T t = (plane.distance - plane.normal.dot(ray.origin)) / denom;
        if(t < ray.tmin || t > ray.tmax) {
            return std::nullopt;
        }

        RayHit<T> hit;
        hit.t = t;
        hit.point = ray.origin + ray.direction * t;
        hit.normal = plane.normal; // unit by Plane contract
        hit.backface = denom > T(0); // crossing toward the negative side
        return hit;
    }

} // namespace cad