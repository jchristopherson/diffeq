# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_base.f90
    ${dir}/diffeq_errors.f90
    ${dir}/diffeq_runge_kutta.f90
    ${dir}/diffeq_dprk45_constants.f90
    ${dir}/diffeq_bsrk32_constants.f90
    ${dir}/diffeq_rk853_constants.f90
    ${dir}/diffeq_rosenbrock_constants.f90
    ${dir}/diffeq_stiffly_stable.f90
)
set(DIFFEQ_SOURCES ${DIFFEQ_SOURCES} PARENT_SCOPE)