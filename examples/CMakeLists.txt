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

# Fixed-Step Exponential Example
add_executable(
    expf_example
    expf_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(expf_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(expf_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Adams Method - Fixed - Example
add_executable(
    adams_fixed_example
    adams_fixed_example.f90
    ${PROJECT_SOURCE_DIR}/test/diffeq_models.f90
)
link_library(adams_fixed_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(adams_fixed_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
