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

# Get Dependencies
find_package(fplot QUIET)

if (fplot_FOUND)
    set(fplot_LIBRARY fplot::fplot)
else()
    include(FetchContent)
    FetchContent_Declare(
        fplot
        GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
    )
    FetchContent_MakeAvailable(fplot)
    set(fplot_LIBRARY fplot)
endif()

# RK45 Example
add_executable(
    rk45_example
    rk45_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rk45_example diffeq ${fplot_LIBRARY})

# RK853 Example
add_executable(
    rk853_example
    rk853_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rk853_example diffeq ${fplot_LIBRARY})

# Comparison Example
add_executable(
    comparison_example
    comparison_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(comparison_example diffeq ${fplot_LIBRARY})

# Rosenbrock Example
add_executable(
    rosenbrock_example
    rosenbrock_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(rosenbrock_example diffeq ${fplot_LIBRARY})

# Parametric Example
add_executable(
    parametric_example
    parametric_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(parametric_example diffeq ${fplot_LIBRARY})

# DAE Example
add_executable(
    dae_example
    dae_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(dae_example diffeq ${fplot_LIBRARY})

# Stiff Integrator Comparison Example
add_executable(
    stiff_integrator_example
    stiff_integrator_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(stiff_integrator_example diffeq ${fplot_LIBRARY})

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET stiff_integrator_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:stiff_integrator_example> 
        $<TARGET_FILE_DIR:stiff_integrator_example>
        COMMAND_EXPAND_LISTS
    )
endif()