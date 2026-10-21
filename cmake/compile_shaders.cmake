# ----------------------------------------------------------------------------
# cad_compile_shaders(<target> SOURCES <glsl files...> OUT_DIR <dir> GEN_DIR <dir>)
#
# Compiles GLSL -> SPIR-V with glslc (ships with the LunarG Vulkan SDK) and
# generates GEN_DIR/shaders_gen.hpp, which embeds every module as a constexpr
# byte array. The executable therefore never loads shader files at runtime.
#
# SOURCES may be absolute or relative to the calling listfile.
# ----------------------------------------------------------------------------
function(cad_compile_shaders target)
    set(oneValue OUT_DIR GEN_DIR)
    set(multiValue SOURCES)
    cmake_parse_arguments(CS "" "${oneValue}" "${multiValue}" ${ARGN})

    find_program(GLSLC_EXECUTABLE
        NAMES glslc
        HINTS "$ENV{VULKAN_SDK}/Bin" "$ENV{VULKAN_SDK}/Bin32"
        DOC "glslc from the LunarG Vulkan SDK")
    if(NOT GLSLC_EXECUTABLE)
        message(FATAL_ERROR
            "glslc was not found. Install the LunarG Vulkan SDK and make sure the "
            "VULKAN_SDK environment variable points to it (restart Visual Studio "
            "after installing).")
    endif()
    message(STATUS "glslc: ${GLSLC_EXECUTABLE}")

    file(MAKE_DIRECTORY "${CS_OUT_DIR}" "${CS_GEN_DIR}")

    set(_spv "")
    foreach(_src IN LISTS CS_SOURCES)
        if(NOT IS_ABSOLUTE "${_src}")
            set(_src "${CMAKE_CURRENT_SOURCE_DIR}/${_src}")
        endif()
        get_filename_component(_name  "${_src}" NAME_WE)     # "mesh"
        get_filename_component(_ext   "${_src}" LAST_EXT)    # ".vert"
        string(SUBSTRING "${_ext}" 1 -1 _stage)              # "vert"
        set(_spv_file "${CS_OUT_DIR}/${_name}_${_stage}.spv")

        add_custom_command(
            OUTPUT  "${_spv_file}"
            COMMAND ${GLSLC_EXECUTABLE} -O --target-env=vulkan1.3
                    -o "${_spv_file}" "${_src}"
            DEPENDS "${_src}"
            COMMENT "glslc ${_src} -> ${_spv_file}"
            VERBATIM)

        list(APPEND _spv "${_spv_file}")
    endforeach()

    set(_gen_hpp "${CS_GEN_DIR}/shaders_gen.hpp")
    add_custom_command(
        OUTPUT  "${_gen_hpp}"
        COMMAND ${CMAKE_COMMAND}
                "-DSOURCES=${_spv}"
                "-DOUT=${_gen_hpp}"
                -P "${CMAKE_SOURCE_DIR}/cmake/EmbedShaders.cmake"
        DEPENDS ${_spv}
        COMMENT "embedding SPIR-V modules into ${_gen_hpp}"
        VERBATIM)

    add_custom_target(${target}_shaders DEPENDS "${_gen_hpp}")
    add_dependencies(${target} ${target}_shaders)
endfunction()