# Get Dependencies
include(FetchContent)
FetchContent_Declare(
    fplot
    GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
)
FetchContent_MakeAvailable(fplot)
set(fplot_LIBRARY fplot)

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

# Stiff Integrator Comparison Example
add_executable(
    stiff_integrator_example
    stiff_integrator_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
target_link_libraries(stiff_integrator_example diffeq ${fplot_LIBRARY})
