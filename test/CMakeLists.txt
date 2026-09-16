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

set(diffeq_test_sources
    diffeq_test.f90
    diffeq_models.f90
    diffeq_jacobian_tests.f90
    diffeq_test_runge_kutta.f90
    diffeq_test_implicit_rk.f90
    diffeq_test_bdf.f90
    diffeq_test_pece.f90
)

find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    set(_linalg_build_testing ${BUILD_TESTING})
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(fortran_test_helper)
    set(BUILD_TESTING ${_linalg_build_testing})
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

add_executable(diffeq_tests ${diffeq_test_sources})
target_link_libraries(
    diffeq_tests
    diffeq
    ${FTH_LIBRARIES}
)
add_test(
    NAME diffeq_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:diffeq_tests>
)

add_executable(diffeq_error_tests diffeq_error_tests.f90)
target_link_libraries(diffeq_error_tests diffeq)

foreach(error_case missing_ode short_grid step_limit tolerance_size)
    add_test(
        NAME diffeq_error_${error_case}
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:diffeq_error_tests> ${error_case}
    )
    set_tests_properties(diffeq_error_${error_case} PROPERTIES WILL_FAIL TRUE)
endforeach()

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET diffeq_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_RUNTIME_DLLS:diffeq_tests>
        $<TARGET_FILE_DIR:diffeq_tests>
        COMMAND_EXPAND_LISTS
    )
endif()

