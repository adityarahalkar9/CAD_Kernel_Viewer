# ============================================================================
# GoogleTest suite for the header-only kernel. Pure CPU - no window, no GPU,
# no Vulkan: runs anywhere, including CI.
# ============================================================================
include(FetchContent)

FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Test files are added incrementally; CONFIGURE_DEPENDS discovers them.
file(GLOB CAD_TEST_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/*Tests.cpp")

if(NOT CAD_TEST_SOURCES)
    message(STATUS "tests: no test files yet - test target skipped")
    return()
endif()

add_executable(cad_kernel_tests ${CAD_TEST_SOURCES})

target_link_libraries(cad_kernel_tests PRIVATE
    cad::kernel
    GTest::gtest_main)

cad_target_warnings(cad_kernel_tests)

include(GoogleTest)
gtest_discover_tests(cad_kernel_tests)