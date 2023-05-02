# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(DIFFEQ_SOURCES
    ${dir}/diffeq.f90
    ${dir}/diffeq_ode_container.f90
    ${dir}/diffeq_fs_integrator.f90
    ${dir}/diffeq_rkfixed.f90
    ${dir}/diffeq_rk4.f90
    ${dir}/diffeq_expfixed.f90
    ${dir}/diffeq_multistep_fixed.f90
    ${dir}/diffeq_abmf.f90
    ${dir}/diffeq_vs_integrator.f90
    ${dir}/diffeq_rkvariable.f90
)
set(DIFFEQ_SOURCES ${DIFFEQ_SOURCES} PARENT_SCOPE)