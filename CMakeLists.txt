# This file is part of diffeq.
# 
# diffeq is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# diffeq is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with diffeq. If not, see <https://www.gnu.org/licenses/>.

# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    diffeq
    LANGUAGES Fortran
    VERSION 1.4.0
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STNADARD_REQUIRED TRUE)

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)

# Install
install(
    EXPORT ${PROJECT_NAME}-targets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build DIFFEQ tests: ${BUILD_TESTING}")
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