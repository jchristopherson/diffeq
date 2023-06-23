# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Fetch the proper content
FetchContent_Declare(
    nonlin
    GIT_REPOSITORY "https://github.com/jchristopherson/nonlin"
)

FetchContent_MakeAvailable(nonlin)

if (WIN32)
    if (BUILD_SHARED_LIBS)
        add_custom_command(
            TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:nonlin>
            $<TARGET_FILE_DIR:${PROJECT_NAME}
        )
    endif()
endif()

set(nonlin_INCLUDE_DIR ${nonlin_BINARY_DIR}/include)
set(nonlin_INCLUDE_DIR ${nonlin_INCLUDE_DIR} PARENT_SCOPE)

# Make a parent-scope variable for the library
set(nonlin_LIBRARY nonlin)
set(nonlin_LIBRARY ${nonlin_LIBRARY} PARENT_SCOPE)
