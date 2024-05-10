# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_base.f90
    ${dir}/diffeq_errors.f90
    ${dir}/diffeq_runge_kutta.f90
    ${dir}/diffeq_dprk45_constants.f90
    # ${dir}/diffeq_fixed_step.f90
    # ${dir}/diffeq_rk_fixed_step.f90
    # ${dir}/diffeq_rk4.f90
    # ${dir}/diffeq_multistep_fixed.f90
    # ${dir}/diffeq_abmf.f90
    # ${dir}/diffeq_variable_step.f90
    # ${dir}/diffeq_variable_singlestep.f90
    # ${dir}/diffeq_runge_kutta.f90
    # ${dir}/diffeq_dprk45_constants.f90
    # ${dir}/diffeq_bsrk32_constants.f90
    # ${dir}/diffeq_implicit_runge_kutta.f90
    # ${dir}/diffeq_sdirk4_constants.f90
)
set(DIFFEQ_SOURCES ${DIFFEQ_SOURCES} PARENT_SCOPE)