# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_ode_container.f90
    ${dir}/diffeq_fs_integrator.f90
    ${dir}/diffeq_rkfixed.f90
    ${dir}/diffeq_rk4.f90
    ${dir}/diffeq_multistep_fixed.f90
    ${dir}/diffeq_abmf.f90
    ${dir}/diffeq_vs_integrator.f90
    ${dir}/diffeq_vssi.f90
    ${dir}/diffeq_rkvariable.f90
    ${dir}/diffeq_dprk45.f90
    ${dir}/diffeq_bsrk32.f90
    ${dir}/diffeq_harmonics.f90
    ${dir}/diffeq_harmonics_routines.f90
    ${dir}/diffeq_implicit_rk.f90
    ${dir}/diffeq_dirk.f90
    ${dir}/diffeq_sdirk4.f90
)
set(DIFFEQ_SOURCES ${DIFFEQ_SOURCES} PARENT_SCOPE)