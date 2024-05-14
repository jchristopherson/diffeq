include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# RK45 Example
add_executable(
    rk45_example
    rk45_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(rk45_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rk45_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# RK853 Example
add_executable(
    rk853_example
    rk853_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(rk853_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rk853_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Comparison Example
add_executable(
    comparison_example
    comparison_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(comparison_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})

# Rosenbrock Example
add_executable(
    rosenbrock_example
    rosenbrock_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(rosenbrock_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rosenbrock_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
