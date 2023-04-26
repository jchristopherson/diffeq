include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)

# RK4 Example
add_executable(
    rk4_example
    rk4_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(rk4_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rk4_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})