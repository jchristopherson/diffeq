!> @brief This module contains several ODE solvers and associated types.
module diffeq
    use iso_fortran_env
    use diffeq_base
    use diffeq_runge_kutta
    use diffeq_stiffly_stable
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
    public :: stiff_single_step_integrator
    public :: form_stiff_single_matrices
    public :: runge_kutta_45
    public :: runge_kutta_23
    public :: runge_kutta_853
    public :: rosenbrock
    
    
    ! public :: DIFFEQ_MEMORY_ALLOCATION_ERROR
    ! public :: DIFFEQ_NULL_POINTER_ERROR
    ! public :: DIFFEQ_MATRIX_SIZE_ERROR
    ! public :: DIFFEQ_ARRAY_SIZE_ERROR
    ! public :: DIFFEQ_INVALID_INPUT_ERROR
    ! public :: DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
    ! public :: DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR
    ! public :: DIFFEQ_INVALID_OPERATION_ERROR

end module