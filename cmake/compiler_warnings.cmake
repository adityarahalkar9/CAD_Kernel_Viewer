# ----------------------------------------------------------------------------
# cad_target_warnings(<target>)
#
# Applies the project's warning profile to a target built by *this* project.
# Third-party dependencies (SDL3, GoogleTest) are deliberately left alone.
# ----------------------------------------------------------------------------
function(cad_target_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE
            /W4                 # high warning level
            /permissive-        # strict standard conformance
            /utf-8              # source & execution charset
            /Zc:__cplusplus     # correct __cplusplus macro
            /MP                 # parallel compilation
            /wd4127             # conditional expression is constant (STL/gtest macros)
            /wd4324)            # structure padding due to alignment (math types)
        if(CAD_ENABLE_WERROR)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow -Wold-style-cast)
        if(CAD_ENABLE_WERROR)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()
endfunction()