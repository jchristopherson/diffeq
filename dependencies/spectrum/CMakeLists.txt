# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    spectrum
    GIT_REPOSITORY "https://github.com/jchristopherson/spectrum"
    GIT_TAG main
)

FetchContent_MakeAvailable(spectrum)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:spectrum>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(spectrum_INCLUDE_DIR ${spectrum_BINARY_DIR}/include)
set(spectrum_INCLUDE_DIR ${spectrum_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(spectrum_LIBRARY spectrum)
set(spectrum_LIBRARY ${spectrum_LIBRARY} PARENT_SCOPE)
