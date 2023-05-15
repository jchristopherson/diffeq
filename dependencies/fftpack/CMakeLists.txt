# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    fftpack
    GIT_REPOSITORY "https://github.com/fortran-lang/fftpack"
    GIT_TAG main
)

FetchContent_MakeAvailable(fftpack)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:fftpack>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(fftpack_INCLUDE_DIR ${fftpack_BINARY_DIR}/include)
set(fftpack_INCLUDE_DIR ${fftpack_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(fftpack_LIBRARY fftpack)
set(fftpack_LIBRARY ${fftpack_LIBRARY} PARENT_SCOPE)
