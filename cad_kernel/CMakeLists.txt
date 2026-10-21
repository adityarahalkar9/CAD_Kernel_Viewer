# ============================================================================
# The kernel is deliberately HEADER-ONLY:
#   * zero build friction - consumers just #include it
#   * usable from tests, tools, and future CAD services alike
#   * templates are instantiated in each translation unit as needed
# ============================================================================

add_library(cad_kernel INTERFACE)
add_library(cad::kernel ALIAS cad_kernel)

target_include_directories(cad_kernel INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_compile_features(cad_kernel INTERFACE cxx_std_20)

# List every header of the final design. Files that do not exist yet (the
# project is assembled incrementally) are silently skipped, and appear in the
# IDE automatically once created and CMake re-runs.
set(CAD_KERNEL_HEADERS
    cad_kernel.hpp
    support/Assertions.hpp
    support/Log.hpp
    math/constants.hpp
    math/utils.hpp
    math/vec2.hpp
    math/vec3.hpp
    math/vec4.hpp
    math/mat3.hpp
    math/mat4.hpp
    math/quaternion.hpp
    math/transform.hpp
    math/geometry_types.hpp
    geometry/ray_intersections.hpp
    geometry/closest_point.hpp
    geometry/plane_ops.hpp
    mesh/handles.hpp
    mesh/half_edge_mesh.hpp
    mesh/mesh_builder.hpp
    mesh/mesh_algorithms.hpp
    mesh/mesh_io.hpp
    spatial/AABBTree.hpp
    spatial/mesh_bvh.hpp
)

set(_present "")
foreach(_h IN LISTS CAD_KERNEL_HEADERS)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_h}")
        list(APPEND _present "${_h}")
    endif()
endforeach()

if(_present)
    target_sources(cad_kernel INTERFACE
        FILE_SET HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        FILES ${_present})
endif()