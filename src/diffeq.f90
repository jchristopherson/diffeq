!> @brief This module contains several ODE solvers and associated types.
module diffeq
    use iso_fortran_env
    use diffeq_base
    use diffeq_runge_kutta
    ! use diffeq_fixed_step
    ! use diffeq_abmf
    ! use diffeq_rk_fixed_step
    ! use diffeq_rk4
    ! use diffeq_multistep_fixed
    ! use diffeq_variable_step
    ! use diffeq_variable_singlestep
    ! use diffeq_runge_kutta
    ! use diffeq_implicit_runge_kutta
    implicit none
    private
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
    public :: single_step_interpolate
    public :: single_step_integrator
    public :: runge_kutta_45
    public :: runge_kutta_23
    public :: runge_kutta_853
    
    ! public :: ode
    ! public :: ode_jacobian
    ! public :: ode_mass_matrix
    ! public :: ode_container
    ! public :: ode_integrator
    ! public :: ode_solver
    ! public :: ode_integer_inquiry
    ! public :: fixed_step_integrator
    ! public :: ode_fixed_step
    ! public :: rk_fixed_integrator
    ! public :: rkf_get_array_parameter
    ! public :: rkf_get_matrix_parameter
    ! public :: rk4_fixed_integrator
    ! public :: fixed_multistep_integrator
    ! public :: adams_fixed_integerator
    ! public :: variable_step_integrator
    ! public :: variable_step_attempt
    ! public :: variable_step_action
    ! public :: variable_step_interpolation
    ! public :: next_step_size_calculator
    ! public :: variable_singlestep_integrator
    ! public :: pi_controller
    ! public :: rkv_get_matrix_parameter
    ! public :: rkv_get_array_parameter
    ! public :: rkv_get_boolean_parameter
    ! public :: rkv_get_integer_parameter
    ! public :: rkv_action
    ! public :: rkv_set_up_interp
    ! public :: rk_variable_integrator
    ! public :: dprk45_integrator
    ! public :: bsrk32_integrator
    ! public :: build_factored_newton_matrix_routine
    ! public :: implicit_rk_variable_integrator
    ! public :: sdirk_integrator
    ! public :: sdirk4_integrator
    ! public :: DIFFEQ_MEMORY_ALLOCATION_ERROR
    ! public :: DIFFEQ_NULL_POINTER_ERROR
    ! public :: DIFFEQ_MATRIX_SIZE_ERROR
    ! public :: DIFFEQ_ARRAY_SIZE_ERROR
    ! public :: DIFFEQ_INVALID_INPUT_ERROR
    ! public :: DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
    ! public :: DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR
    ! public :: DIFFEQ_INVALID_OPERATION_ERROR

end module