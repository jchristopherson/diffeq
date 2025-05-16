!> @brief This module contains several ODE solvers and associated types.
module diffeq
    use iso_fortran_env
    use diffeq_base
    use diffeq_runge_kutta
    use diffeq_implicit_runge_kutta
    use diffeq_vode
    implicit none
    private
    
    ! DIFFEQ_BASE.F90
    public :: ode
    public :: ode_jacobian
    public :: ode_mass_matrix
    public :: ode_container
    public :: ode_integrator
    public :: ode_solver
    public :: ode_integer_inquiry
    public :: attempt_single_step
    public :: get_single_step_logical_parameter
    public :: single_step_post_step_routine
    public :: single_step_pre_step_routine
    public :: single_step_interpolate
    public :: single_step_integrator

    ! DIFFEQ_RUNGE_KUTTA.F90
    public :: runge_kutta_45
    public :: runge_kutta_23
    public :: runge_kutta_853

    ! DIFFEQ_IMPLICIT_RUNGE_KUTTA.F90
    public :: rosenbrock

    ! DIFFEQ_VODE.F90
    public :: VODE_ADAMS_METHOD
    public :: VODE_BDF_METHOD
    public :: vode
    public :: adams
    public :: bdf
end module