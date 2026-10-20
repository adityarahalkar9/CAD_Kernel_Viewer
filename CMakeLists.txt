cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project(CadKernelViewer
    VERSION 0.1.0
    DESCRIPTION "Header-only CAD geometry kernel and Vulkan 1.4 / SDL3 desktop viewer"
    LANGUAGES CXX)

# ============================================================================
# Options
# ============================================================================
option(CAD_ENABLE_TESTS      "Build the GoogleTest unit-test suite"              ON)
option(CAD_ENABLE_VALIDATION "Request Vulkan validation layers in Debug builds" ON)
option(CAD_ENABLE_WERROR     "Treat compiler warnings as errors"                 OFF)
option(CAD_USE_SYSTEM_SDL3   "Use an installed SDL3 (vcpkg/system) instead of fetching" OFF)

set(CAD_SDL3_GIT_TAG "release-3.2.0" CACHE STRING
    "SDL3 git tag fetched from GitHub when no system SDL3 is found")

# ============================================================================
# Global language settings
# ============================================================================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Single-config generators (Ninja/make) need a default build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# Third-party projects (SDL3, GoogleTest) are built statically and linked into
# our binaries: the final CadViewer.exe is fully self-contained, no DLLs.
if(NOT DEFINED BUILD_SHARED_LIBS)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompilerWarnings)

# ============================================================================
# Vulkan (requires the LunarG Vulkan SDK, 1.4 or newer)
# ============================================================================
find_package(Vulkan REQUIRED)
message(STATUS "Vulkan SDK: ${Vulkan_VERSION} (${Vulkan_INCLUDE_DIRS})")
if(Vulkan_VERSION VERSION_LESS 1.4)
    message(WARNING
        "The installed Vulkan SDK is ${Vulkan_VERSION}, but this project targets "
        "Vulkan 1.4. Install the latest LunarG SDK, then restart Visual Studio so "
        "the VULKAN_SDK environment variable is visible to CMake.")
endif()

# ============================================================================
# Subprojects
#   The project is assembled incrementally (kernel first, then viewer, then
#   tests). Each subdirectory is only added once its CMakeLists.txt exists,
#   and each of them tolerates having no sources yet.
# ============================================================================
add_subdirectory(cad_kernel)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/viewer/CMakeLists.txt")
    add_subdirectory(viewer)
else()
    message(STATUS "viewer/ not present yet - skipping viewer application")
endif()

if(CAD_ENABLE_TESTS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
    enable_testing()
    add_subdirectory(tests)
endif()