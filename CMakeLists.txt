# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    diffeq
    LANGUAGES Fortran
    VERSION 1.1.1
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
add_subdirectory(dependencies)

# Source
add_subdirectory(src)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${DIFFEQ_SOURCES}
)
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${linalg_LIBRARY} ${linalg_INCLUDE_DIR})

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Examples
option(BUILD_DIFFEQ_EXAMPLES "Build DIFFEQ examples")
message(STATUS "Build DIFFEQ examples: ${BUILD_DIFFEQ_EXAMPLES}")
if (BUILD_DIFFEQ_EXAMPLES)
    add_subdirectory(examples)
endif()