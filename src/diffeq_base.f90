! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module diffeq_base
    !! A collection of base types for the DIFFEQ library.
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use diffeq_errors
    use diffeq_containers, only : ode_results
    implicit none
    private
    public :: ode
    public :: ode_jacobian
    public :: ode_mass_matrix
    public :: ode_container
    public :: ode_integrator
    public :: ode_solver
    public :: ode_integer_inquiry
    public :: ode_dense_output
    public :: ode_event_function
    public :: ode_event_solver
    public :: ode_event_statistics
    
! ------------------------------------------------------------------------------
    interface
        subroutine ode(x, y, dydx, args)
            !! An interface to a routine containing a system of first order
            !! ODE's to integrate.
            use iso_fortran_env
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the
                !! dependent variables.
            real(real64), intent(out), dimension(:) :: dydx
                !! An N-element array where the derivative values will be
                !! written.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the routine.
        end subroutine

        subroutine ode_jacobian(x, y, jac, args)
            !! An interface to a routine capable of computing the Jacobian
            !! matrix of a system of first order ODE's.
            use iso_fortran_env
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the
                !! dependent variables.
            real(real64), intent(out), dimension(:,:) :: jac
                !! An M-by-N matrix where the Jacobian will be written.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the routine.
        end subroutine

        subroutine ode_mass_matrix(x, y, m, args)
            !! An interface to a routine capable of computing a mass matrix
            !! for a system of differential equations.
            use iso_fortran_env
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the
                !! dependent variables.
            real(real64), intent(out), dimension(:,:) :: m
                !! An N-by-N matrix where the 
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the routine.
        end subroutine

        subroutine ode_event_function(x, y, f, args)
            !! A callback routine for evaluating event status.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current state vector.
            real(real64), intent(out), dimension(:) :: f
                !! An M-element array containing the values of the event
                !! functions.  The solver attempts to find the situation 
                !! where f equals zero thereby signaling an event.  Notice, M 
                !! need not equal N (N is the number of state variables).
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the routine.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type ode_container
        !! Defines the callbacks and model options for an ODE problem.
        !!
        !! The required right-hand-side callback describes
        !! \( y' = f(x,y) \). Optional callbacks provide an analytical
        !! Jacobian \( J = \partial f / \partial y \) and a mass matrix for
        !! systems of the form \( M(x,y)y' = f(x,y) \). If no Jacobian is
        !! supplied, DIFFEQ estimates it using finite differences. A null
        !! mass-matrix callback means that \( M = I \).
        logical, private :: m_massDependent = .false.
            ! A value determining if the mass matrix is state dependent such 
            ! that it must be recomputed at each step. 
        real(real64), private, allocatable, dimension(:) :: m_jwork
            ! Jacobian calculation workspace array.
        real(real64), private :: m_fdStep = sqrt(epsilon(1.0d0))
            ! Finite difference step size.
        procedure(ode), pointer, public, nopass :: fcn => null()
            !! A pointer to the routine containing the ODEs to integrate.
        procedure(ode_jacobian), pointer, public, nopass :: &
            jacobian => null()
            !! A pointer to the routine containing the analytical Jacobian.
            !! If supplied, this routine is utilized; however, if null, a finite
            !! difference approximation is utilized.
        procedure(ode_mass_matrix), pointer, public, nopass :: &
            mass_matrix => null()
            !! A pointer to the routine containing the mass matrix for the
            !! system.  If set to null (the default), an identity mass matrix 
            !! will be assumed.
        procedure(ode_event_function), pointer, public, nopass :: &
            event_function => null()
            !! A pointer to the event function.
        integer, private :: m_nevents = 0
            !! The number of event equations.  This number need not be equal
            !! to the number of ODE's.
    contains
        procedure, private :: allocate_workspace => oc_alloc_workspace
        procedure, public :: get_finite_difference_step => oc_get_fd_step
        procedure, public :: set_finite_difference_step => oc_set_fd_step
        procedure, public :: get_is_mass_matrix_dependent => &
            oc_get_is_mass_dependent
        procedure, public :: set_is_mass_matrix_dependent => &
            oc_set_is_mass_dependent
        procedure, public :: compute_jacobian => oc_jacobian
        procedure, public :: get_is_ode_defined => oc_get_is_ode_defined
        procedure, public :: use_finite_difference_jacobian => &
            oc_use_finite_difference_jacobian
        procedure, public :: compute_rhs => oc_eval_rhs
        procedure, public :: compute_mass_matrix => oc_eval_mass_matrix
        procedure, public :: has_mass_matrix => oc_has_mass_matrix
        procedure, public :: get_is_ode_event_function_defined => &
            oc_is_event_fcn_defined
        procedure, public :: get_event_equation_count => oc_get_event_eqn_count
        procedure, public :: set_event_equation_count => oc_set_event_eqn_count
    end type

! ------------------------------------------------------------------------------
    type ode_event_statistics
        !! An ODE event solver statistics container.
        integer :: interation_count = 0
            !! The number of iterations taken by the solver.
        integer :: function_evaluations = 0
            !! The number of event function evaluations.
        integer :: event_index = 0
            !! The index of the event equation where an event was located.
    end type

    type ode_event_solver
        !! Provides a mechanism for event detection.
        integer, private :: m_maxIter = 1000
            !! The maximum number of iterations allowed.
        real(real64), private :: m_tol = 1.0d-12
            !! The convergence tolerance.
        type(ode_event_statistics), public :: solver_statistics
            !! Event solver information.
    contains
        procedure, public :: get_max_allowed_iterations => oes_get_max_iter
            !! Gets the number of iterations the solver is allowed to take.
        procedure, public :: set_max_allowed_iterations => oes_set_max_iter
            !! Sets the number of iterations the solver is allowed to take.
        procedure, public :: get_solver_tolerance => oes_get_tol
            !! Gets the solver tolerance.
        procedure, public :: set_solver_tolerance => oes_set_tol
            !! Sets the solver tolerance.
        procedure, public :: is_event_suspected => oes_event_in_bracket
            !! Determines if an event is suspected to exist within the &
            !! specified bracket.
        procedure, public :: find_event => oes_find_event
            !! Attempts to find any events within a specified bracket.
        procedure, private :: find_channel_event => oes_find_single_channel_event
    end type

! ------------------------------------------------------------------------------
    type, abstract :: ode_integrator
        !! Abstract base class for all DIFFEQ ODE integrators.
        !!
        !! Every concrete integrator accepts an `ode_container`, an array of
        !! output coordinates, and an initial state. Adaptive solvers use a
        !! scaled local error such as
        !! \( \mathrm{err}_i / (\mathrm{atol} + \mathrm{rtol}|y_i|) \)
        !! to choose their step sizes.
        type(ode_results), private :: m_buffer
            ! The internal solution buffer.
        real(real64), private :: m_abstol = 1.0d-6
            ! The absolute tolerance value applied to each equation.
        real(real64), private :: m_reltol = 1.0d-6
            ! The relative tolerance value applied to each equation.
        real(real64), private, allocatable, dimension(:) :: m_abstols
            ! Per-equation absolute tolerances, when defined.
        real(real64), private, allocatable, dimension(:) :: m_reltols
            ! Per-equation relative tolerances, when defined.
        real(real64), private :: m_minStep = 1.0d1 * epsilon(1.0d0)
            ! The minimum allowable step size.
        real(real64), private :: m_maxStep = huge(1.0d0)
            ! The maximum allowable step size.
        real(real64), private :: m_safetyFactor = 0.9d0
            ! The step size safety factor.
        real(real64), private :: m_beta = 0.0d0
            ! PI step size controller exponent.
        logical, private :: m_reject = .false.
            ! Internal variable tracking step size acceptance.
        integer, private :: m_stepLimit = 1000000
            ! A limit on the total number of integration steps.  Exceeding this
            ! value may indicate that a different integrator should be
            ! considered.  Either that, or this is a really large-sized
            ! problem that is being solved.
        logical, private :: m_allowOvershoot = .true.
            ! True if the solver is allowed to overshoot the final value and
            ! interpolate back.  False, if the solver must terminate on the
            ! final value.
        logical, private :: m_denseEnabled = .false.
            ! Dense output enabled?
        logical, private :: m_denseValid = .false.
            ! Is dense output valid for this solver type?
        real(real64), private :: m_denseInterval(2) = 0.0d0
            ! The interval over which the stored dense output interpolation
            ! is defined
        integer, private :: m_denseSize = 0
            ! The number of dense output equations.
        type(ode_event_solver), public :: event_solver
            !! The ODE event solver.  Enable event detection by defining the
            !! event_function in the [[ode_controller]] object.
    contains
        procedure, public :: set_dense_output_enabled => oi_set_dense_enabled
            !! Enables dense output on every accepted step (default: false).
        procedure, public :: get_dense_output_enabled => oi_get_dense_enabled
            !! Gets the explicit dense-output preparation setting.
        procedure, public :: has_dense_output => oi_has_dense_output
            !! Tests whether a prepared accepted-step interpolant is available.
        procedure, public :: get_dense_interval => oi_get_dense_interval
            !! Gets the prepared step's start and end in integration order.
        procedure, public :: evaluate_dense => oi_evaluate_dense
            !! Evaluates the prepared step without advancing integration.
        procedure(ode_dense_output), public, pass, deferred :: compute_dense_output
            !! Implementation hook for evaluating a validated dense-output query.
        procedure, public :: set_dense_interval => oi_set_dense_interval
            !! Implementation hook publishing a prepared accepted step.
        procedure, public :: invalidate_dense_output => oi_invalidate_dense
            !! Invalidates the prepared interpolant before stepping or resetting.
        procedure(ode_solver), public, pass, deferred :: solve
            !! Solves the supplied system of ODE's.
        procedure(ode_integer_inquiry), public, pass, deferred :: get_order
            !! Returns the order of the integrator.
        procedure, public :: append_to_buffer => oi_append_to_buffer
            !! Appends the supplied solution point to the internal solution 
            !! buffer.
        procedure, public :: get_solution => oi_get_solution
            !! Returns the solution computed by the integrator.
        procedure, public :: get_ode_results => oi_get_ode_results
            !! Returns the solution stored within an [[ode_results]] object.
        procedure, public :: clear_buffer => oi_clear_buffer
            !! Clears the contents of the buffer.
        procedure, public :: get_absolute_tolerance => oi_get_abs_tol
            !! Gets the absolute error tolerance.
        procedure, private :: oi_set_abs_tol
        procedure, private :: oi_set_abs_tols
        generic, public :: set_absolute_tolerance => oi_set_abs_tol, &
            oi_set_abs_tols
            !! Sets the absolute error tolerance.
        procedure, public :: get_relative_tolerance => oi_get_rel_tol
            !! Gets the relative error tolerance.
        procedure, private :: oi_set_rel_tol
        procedure, private :: oi_set_rel_tols
        generic, public :: set_relative_tolerance => oi_set_rel_tol, &
            oi_set_rel_tols
            !! Sets the relative error tolerance.
        procedure, public :: compute_error_norm => oi_estimate_error
            !! Computes the norm of the scaled error estimate.
        procedure, public :: get_minimum_step_size => oi_get_min_step
            !! Gets the magnitude of the minimum allowed step size.
        procedure, public :: set_minimum_step_size => oi_set_min_step
            !! Sets the magnitude of the minimum allowed step size.
        procedure, public :: get_maximum_step_size => oi_get_max_step
            !! Gets the magnitude of the maximum allowed step size.
        procedure, public :: set_maximum_step_size => oi_set_max_step
            !! Sets the magnitude of the maximum allowed step size.
        procedure, public :: get_step_size_factor => oi_get_safety_factor
            !! Gets the step size safety factor.
        procedure, public :: set_step_size_factor => oi_set_safety_factor
            !! Sets the step size safety factor.
        procedure, public :: get_step_size_control_parameter => &
            oi_get_control_parameter
            !! Gets the step size PI control parameter.
        procedure, public :: set_step_size_control_parameter => &
            oi_set_control_parameter
            !! Sets the step size PI control parameter.
        procedure, public :: estimate_next_step_size => oi_next_step
            !! Estimates the next step size.
        procedure, public :: estimate_inital_step_size => oi_initial_step
            !! Computes an estimate of an initial step size.
        procedure, public :: get_step_limit => oi_get_step_limit
            !! Gets the limit on the number of integration steps.
        procedure, public :: set_step_limit => oi_set_step_limit
            !! Sets the limit on the number of integration steps.
        procedure, public :: get_allow_overshoot => oi_get_allow_overshoot
            !! Gets a value determining if the solver is allowed to overshoot 
            !! the final value in the integration range.
        procedure, public :: set_allow_overshoot => oi_set_allow_overshoot
            !! Sets a value determining if the solver is allowed to overshoot 
            !! the final value in the integration range.
        procedure, public :: log_failed_step => oi_log_failed_step
            !! Logs a failed integration step.
        procedure, public :: log_successful_step => oi_log_successful_step
            !! Logs a successful integration step.
        procedure, public :: log_jacobian_evaluation => oi_log_jacobian_eval
            !! Logs a Jacobian matrix evaluation.
        procedure, public :: log_function_evaluation => oi_log_function_eval
            !! Logs a function evaluation.
        procedure, public :: log_matrix_factorization => oi_log_matrix_factorization
            !! Logs a matrix factorization.
    end type

    abstract interface
        subroutine ode_dense_output(this, x, y)
            import ode_integrator, real64
            class(ode_integrator), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64), intent(out) :: y(:)
        end subroutine

        subroutine ode_solver(this, sys, x, iv, args)
            !! Solves the supplied system of ODE's.
            use iso_fortran_env
            import ode_integrator
            import ode_container
            class(ode_integrator), intent(inout) :: this
                !! The ode_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in), dimension(:) :: x
                !! An array, of at least 2 values, defining at a minimum
                !! the starting and ending values of the independent variable 
                !! integration range.  If more than two values are specified, 
                !! the integration results will be returned at the supplied 
                !! values.
            real(real64), intent(in), dimension(:) :: iv
                !! An array containing the initial values for each ODE.
            class(*), intent(inout), optional, target :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        pure function ode_integer_inquiry(this) result(rst)
            !! Returns an integer value from the ode_integrator object.
            use iso_fortran_env
            import ode_integrator
            class(ode_integrator), intent(in) :: this
                !! The ode_integrator object.
            integer :: rst
                !! The requested value.
        end function
    end interface

contains


! ******************************************************************************
! ODE_CONTAINER ROUTINES
! ------------------------------------------------------------------------------
pure function oc_get_fd_step(this) result(rst)
    !! Gets the size of the step to use for the finite difference
    !! calculations used to estimate the Jacobian.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    real(real64) :: rst
        !! The step size.
    rst = this%m_fdStep
end function

! --------------------
subroutine oc_set_fd_step(this, x)
    !! Sets the size of the step to use for the finite difference
    !! calculations used to estimate the Jacobian.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    real(real64), intent(in) :: x
        !! The step size.
    this%m_fdStep = x
end subroutine

! ------------------------------------------------------------------------------
pure function oc_get_is_mass_dependent(this) result(rst)
    !! Gets a value determining if the mass matrix is state-dependent
    !! such that it requires updating at every integration step.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    logical :: rst
        !! True if the mass matrix is state-dependent such that it 
        !! requires updating at each integration step; else, false if the
        !! mass matrix is not state-dependent and can be treated as constant
        !! for all integration steps.
    rst = this%m_massDependent
end function

! --------------------
subroutine oc_set_is_mass_dependent(this, x)
    !! Sets a value determining if the mass matrix is state-dependent
    !! such that it requires updating at every integration step.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    logical :: x
        !! True if the mass matrix is state-dependent such that it 
        !! requires updating at each integration step; else, false if the
        !! mass matrix is not state-dependent and can be treated as constant
        !! for all integration steps.
    this%m_massDependent = x
end subroutine

! ------------------------------------------------------------------------------
subroutine oc_jacobian(this, x, y, jac, args)
    !! Computes the Jacobian matrix for the system of ODEs.  If
    !! a routine is provided with an analytical Jacobian, the supplied
    !! routine is utilized; else, the Jacobian is estimated via a forward
    !! difference approximation.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    real(real64), intent(in) :: x
        !! The current independent variable value.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current dependent
        !! variable values.
    real(real64), intent(out), dimension(:,:) :: jac
        !! An N-by-N matrix where the Jacobian will be written.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: i, ndof
    real(real64) :: h
    
    ! Initialization
    ndof = size(y)
    h = this%get_finite_difference_step()

    ! Input Checking
    if (size(jac, 1) /= ndof .or. size(jac, 2) /= ndof) &
        error stop DIFFEQ_MATRIX_SIZE_ERROR

    ! Use a user-defined routine, and then be done
    if (associated(this%jacobian)) then
        call this%jacobian(x, y, jac, args)
        return
    end if

    ! Allocate workspace.  No action is taken if the proper workspace is
    ! already allocated.
    call this%allocate_workspace(ndof)

    ! Finite Difference Approximation
    ! J(i,j) = df(i) / dy(j)
    this%m_jwork(1:ndof) = y
    call this%compute_rhs(x, y, this%m_jwork(ndof+1:), args)
    do i = 1, ndof
        this%m_jwork(i) = this%m_jwork(i) + h
        call this%compute_rhs(x, this%m_jwork(1:ndof), jac(:,i), args)
        jac(:,i) = (jac(:,i) - this%m_jwork(ndof+1:)) / h
        this%m_jwork(i) = y(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
function oc_use_finite_difference_jacobian(this) result(rst)
    !! States if a finite difference approximation to the Jacobian is utilized.
    class(ode_container), intent(in) :: this
        !! The [[ode_container]] object.
    logical :: rst
        !! True if a finite difference approximation is being utilized; else,
        !! false.

    rst = .not.associated(this%jacobian)
end function

! ------------------------------------------------------------------------------
subroutine oc_alloc_workspace(this, ndof)
    ! Use to allocate internal workspaces.  This routine only takes action
    ! if the workspace array(s) are not sized properly for the application.
    class(ode_container), intent(inout) :: this
        ! The ode_container object.
    integer, intent(in) :: ndof
        ! The number of degrees of freedom.

    ! Jacobian Workspace Allocation
    if (allocated(this%m_jwork)) then
        if (size(this%m_jwork) /= 2 * ndof) then
            deallocate(this%m_jwork)
            allocate(this%m_jwork(2 * ndof), source = 0.0d0)
        end if
    else
        allocate(this%m_jwork(2 * ndof), source = 0.0d0)
    end if
end subroutine

! ------------------------------------------------------------------------------
function oc_get_is_ode_defined(this) result(rst)
    !! Gets a logical value determining if the ODE routine has been defined.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    logical :: rst
        !! True if the ODE routine has been defined; else, false.

    rst = associated(this%fcn)
end function

! ------------------------------------------------------------------------------
subroutine oc_eval_rhs(this, x, y, dydx, args)
    !! Evaluates the right-hand-side routine defined by the user-defined
    !! function pointer 'fcn'.  If no routine has been assigned to the 
    !! function pointer a DIFFEQ_NULL_POINTER_ERROR is thrown.
    class(ode_container), intent(inout) :: this
        !! The [[ode_container]] object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the
        !! dependent variables.
    real(real64), intent(out), dimension(:) :: dydx
        !! An N-element array where the derivative values will be
        !! written.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the routine.

    if (.not.this%get_is_ode_defined()) error stop DIFFEQ_NULL_POINTER_ERROR
    call this%fcn(x, y, dydx, args = args)
end subroutine

! ------------------------------------------------------------------------------
function oc_has_mass_matrix(this) result(rst)
    !! Returns true if the problem utilizes a mass matrix.
    class(ode_container), intent(in) :: this
        !! The [[ode_container]] object.
    logical :: rst
        !! True if a mass matrix is utilized; else, false.

    rst = associated(this%mass_matrix)
end function

! ------------------------------------------------------------------------------
subroutine oc_eval_mass_matrix(this, x, y, m, args)
    use linalg, only : identity
    !! Evaluates the mass matrix defined by the function pointer 'mass_matrix'.
    !! If no function is assigned to 'mass_matrix' this routine returns an
    !! appropriately sized identity matrix.
    class(ode_container), intent(inout) :: this
        !! The [[ode_container]] object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the
        !! dependent variables.
    real(real64), intent(out), dimension(:,:) :: m
        !! An N-by-N matrix where the 
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the routine.

    integer :: n
    n = size(y)
    if (associated(this%mass_matrix)) then
        call this%mass_matrix(x, y, m, args = args)
    else
        m = identity(n)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function oc_is_event_fcn_defined(this) result(rst)
    !! Tests to see if the event function has been defined.
    class(ode_container), intent(in) :: this
        !! The [[ode_container]] object.
    logical :: rst
        !! True if the event function has been defined; else, false.

    rst = associated(this%event_function)
end function

! ------------------------------------------------------------------------------
pure function oc_get_event_eqn_count(this) result(rst)
    !! Gets the number of event functions.
    class(ode_container), intent(in) :: this
        !! The [[ode_container]] object.
    integer :: rst
        !! The equation count.
    rst = this%m_nevents
end function

! --------------------
subroutine oc_set_event_eqn_count(this, x)
    !! Sets the number of event functions.
    class(ode_container), intent(inout) :: this
        !! The [[ode_container]] object.
    integer, intent(in) :: x
        !! The equation count.
    this%m_nevents = x
end subroutine

! ******************************************************************************
! ODE_INTEGRATOR ROUTINES
! ------------------------------------------------------------------------------
subroutine oi_append_to_buffer(this, x, y)
    !! Appends the supplied solution point to the internal solution buffer.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The independent variable value.
    real(real64), intent(in), dimension(:) :: y
        !! The values of the dependent variables corresponding to x.

    call this%m_buffer%push(x, y)
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_solution(this) result(rst)
    !! Returns the solution computed by the integrator stored as a matrix with
    !! the first column containing the values of the independent variable at
    !! which the solution was computed.  The remaining columns contain the
    !! solutions for each of the integrated equations in the order in which they
    !! appear in the source routine.  Notice, the solve routine must be called
    !! before this routine.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The resulting solution matrix.

    if (this%m_buffer%count() == 0) then
        allocate(rst(0, 0))
    else
        rst = this%m_buffer%to_matrix()
    end if
end function

! ------------------------------------------------------------------------------
pure function oi_get_ode_results(this) result(rst)
    !! Returns a copy of the solution computed as a [[ode_results]] container.
    !! The integrator retains the original copy.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    type(ode_results) :: rst
        !! The [[ode_results]] container.

    rst = this%m_buffer
end function

! ------------------------------------------------------------------------------
subroutine oi_clear_buffer(this)
    !! Clears the contents of the buffer.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.

    call this%invalidate_dense_output()
    call this%m_buffer%clear_all()
    this%m_buffer%statistics%failed_steps = 0
    this%m_buffer%statistics%successful_steps = 0
    this%m_buffer%statistics%jacobian_evaluations = 0
    this%m_buffer%statistics%function_evaluations = 0
    this%m_buffer%statistics%matrix_factorizations = 0
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_abs_tol(this, equation) result(rst)
    !! Gets the absolute error tolerance.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    integer, intent(in), optional :: equation
        !! The equation index. If omitted, the scalar tolerance is returned.
    real(real64) :: rst
        !! The tolerance value.
    if (present(equation) .and. allocated(this%m_abstols)) then
        rst = this%m_abstols(equation)
    else
        rst = this%m_abstol
    end if
end function

! --------------------
subroutine oi_set_abs_tol(this, x)
    !! Sets the absolute error tolerance.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance value.
    this%m_abstol = x
    if (allocated(this%m_abstols)) deallocate(this%m_abstols)
end subroutine

! --------------------
subroutine oi_set_abs_tols(this, x)
    !! Sets the absolute error tolerance for each equation.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in), dimension(:) :: x
        !! The tolerance values.
    this%m_abstols = x
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_rel_tol(this, equation) result(rst)
    !! Gets the relative error tolerance.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    integer, intent(in), optional :: equation
        !! The equation index. If omitted, the scalar tolerance is returned.
    real(real64) :: rst
        !! The tolerance value.
    if (present(equation) .and. allocated(this%m_reltols)) then
        rst = this%m_reltols(equation)
    else
        rst = this%m_reltol
    end if
end function

! --------------------
subroutine oi_set_rel_tol(this, x)
    !! Sets the relative error tolerance.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance value.
    this%m_reltol = x
    if (allocated(this%m_reltols)) deallocate(this%m_reltols)
end subroutine

! --------------------
subroutine oi_set_rel_tols(this, x)
    !! Sets the relative error tolerance for each equation.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in), dimension(:) :: x
        !! The tolerance values.
    this%m_reltols = x
end subroutine

! ------------------------------------------------------------------------------
pure function oi_estimate_error(this, y, yest, yerr) result(rst)
    !! Computes the norm of the scaled error estimate.  A value less than one
    !! indicates a successful step.  A value greater than one suggests that the
    !! results do not meet the requested tolerances.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64), intent(in), dimension(:) :: y
        !! The previously accepted solution array (N-element).
    real(real64), intent(in), dimension(size(y)) :: yest
        !! An N-element array containing the next solution point estimate.
    real(real64), intent(in), dimension(size(y)) :: yerr
        !! An N-element array containing the estimate of error for each
        !! equation.
    real(real64) :: rst
        !! The norm of the scaled error.

    ! Local Variables
    integer :: i, n
    real(real64) :: sf, atol, rtol

    ! Initialization
    n = size(y)
    if (allocated(this%m_abstols)) then
        if (size(this%m_abstols) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR
    end if
    if (allocated(this%m_reltols)) then
        if (size(this%m_reltols) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR
    end if

    ! Process
    rst = 0.0d0
    do i = 1, n
        atol = this%get_absolute_tolerance(i)
        rtol = this%get_relative_tolerance(i)
        sf = atol + rtol * max(abs(y(i)), abs(yest(i)))
        rst = rst + (yerr(i) / sf)**2
    end do
    rst = sqrt(rst / n)
end function

! ------------------------------------------------------------------------------
pure function oi_get_min_step(this) result(rst)
    !! Gets the magnitude of the minimum allowed step size.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The step size limit.
    rst = this%m_minStep
end function

! --------------------
subroutine oi_set_min_step(this, x)
    !! Sets the magnitude of the minimum allowed step size.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The step size limit.
    this%m_minStep = abs(x)
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_max_step(this) result(rst)
    !! Gets the magnitude of the maximum allowed step size.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The step size limit.
    rst = this%m_maxStep
end function

! --------------------
subroutine oi_set_max_step(this, x)
    !! Sets the magnitude of the maximum allowed step size.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The step size limit.
    this%m_maxStep = abs(x)
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_safety_factor(this) result(rst)
    !! Gets the safety factor (step size multiplier) used to provide a measure
    !! of control to the step size estimate such that \( h_{new} = f_{s} h \),
    !! where \( f_{s} \) is this safety factor.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The safety factor.
    rst = this%m_safetyFactor
end function

! --------------------
subroutine oi_set_safety_factor(this, x)
    !! Sets the safety factor (step size multiplier) used to provide a measure
    !! of control to the step size estimate such that \( h_{new} = f_{s} h \),
    !! where \( f_{s} \) is this safety factor.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The safety factor.
    this%m_safetyFactor = x
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_control_parameter(this) result(rst)
    !! Gets the step size control parameter \( \beta \) used for PI control of 
    !! the step size.  A value of 0 provides a default step size controller
    !! (non-PI); however, a nonzero value of \( \beta \) provides PI control 
    !! that improves stability, but comes with a potential for efficiency loss.
    !! A good estimate for a starting point for this parameter is \( \beta =
    !! \frac{0.4}{k} \) where \( k \) is the order of the integrator.
    !!
    !! The PI controller for step size is defined as follows.
    !! $$ h_{n+1} = f_{s} h_{n} e_{n}^{-\alpha} e_{n-1}^{\beta} $$
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The control parameter.
    rst = this%m_beta
end function

! --------------------
subroutine oi_set_control_parameter(this, x)
    !! Sets the step size control parameter \( \beta \) used for PI control of 
    !! the step size.  A value of 0 provides a default step size controller
    !! (non-PI); however, a nonzero value of \( \beta \) provides PI control 
    !! that improves stability, but comes with a potential for efficiency loss.
    !! A good estimate for a starting point for this parameter is \( \beta =
    !! \frac{0.4}{k} \) where \( k \) is the order of the integrator.
    !!
    !! The PI controller for step size is defined as follows.
    !! $$ h_{n+1} = f_{s} h_{n} e_{n}^{-\alpha} e_{n-1}^{\beta} $$
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The control parameter
    this%m_beta = x
end subroutine

! ------------------------------------------------------------------------------
function oi_next_step(this, e, eold, h, x) result(rst)
    !! Estimates the next step size based upon the current and previous error
    !! estimates.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: e
        !! The norm of the current scaled error estimate.
    real(real64), intent(inout) :: eold
        !! The norm of the previous step's scaled error estimate.  On output,
        !! this variable is updated.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current independent variable value.
    real(real64) :: rst
        !! The new step size estimate.

    ! Parameters
    real(real64), parameter :: minscale = 2.0d-1
    real(real64), parameter :: maxscale = 1.0d1

    ! Local Variables
    integer :: k
    real(real64) :: alpha, beta, fs, maxstep, minstep, scale
    
    ! Initialization
    k = this%get_order()
    beta = this%get_step_size_control_parameter()
    alpha = 1.0d0 / k - 0.75d0 * beta
    maxstep = this%get_maximum_step_size()
    minstep = this%get_minimum_step_size()
    fs = this%get_step_size_factor()

    ! Process
    if (e <= 1.0d0) then
        ! The step met error tolerances and is acceptable
        if (e == 0.0d0) then
            scale = maxscale
        else
            scale = fs * e**(-alpha) * eold**beta
            if (scale < minscale) scale = minscale
            if (scale > maxscale) scale = maxscale
        end if
        if (this%m_reject) then
            ! Don't let the step size increase if the last step was rejected
            rst = h * min(scale, 1.0d0)
        else
            rst = h * scale
        end if
        eold = max(e, 1.0d-4)
        this%m_reject = .false.
    else
        ! The step is rejected, reduce the step size
        scale = fs * e**(-alpha)
        scale = max(scale, minscale)
        this%m_reject = .true.
        rst = h * scale
    end if

    ! Check the step size against limits
    if (abs(rst) > maxstep) then
        rst = sign(maxstep, h)
    end if
    if (abs(rst) < minstep) error stop DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
end function

! ------------------------------------------------------------------------------
subroutine oi_initial_step(this, sys, xo, xf, yo, fo, h, args)
    !! Computes an estimate of an initial step size.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: xo
        !! The initial value of the independent variable.
    real(real64), intent(in) :: xf
        !! The final value of the independent variable.
    real(real64), intent(in), dimension(:) :: yo
        !! The initial values of the dependent variables (N-element).
    real(real64), intent(out), dimension(size(yo)) :: fo
        !! An N-element array where the function values at xo will be written.
    real(real64), intent(out) :: h
        !! The initial step size estimate.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: i, n
    real(real64) :: dx, derivative_norm
    real(real64), dimension(size(yo)) :: tol

    ! Use a very basic estimate of an initial step size.  The catch is that a
    ! single function evaluation must be made; however, this is likely needed
    ! in the first place, so no real extra work is necessary.
    n = size(yo)
    if (allocated(this%m_abstols)) then
        if (size(this%m_abstols) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR
    end if
    if (allocated(this%m_reltols)) then
        if (size(this%m_reltols) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR
    end if
    dx = 0.1d0 * (xf - xo)
    do i = 1, n
        tol(i) = max(this%get_absolute_tolerance(i), &
            this%get_relative_tolerance(i))
    end do
    call sys%compute_rhs(xo, yo, fo, args)
    derivative_norm = norm2(fo / tol)
    if (derivative_norm > 0.0d0) then
        h = min(abs(dx), 2.0d0 / derivative_norm)
    else
        h = abs(dx)
    end if
    h = min(h, this%get_maximum_step_size())
    h = sign(h, dx)
    call this%log_function_evaluation()
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_step_limit(this) result(rst)
    !! Gets the limit on the number of integration steps that may be taken 
    !! before the solver terminates.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    integer :: rst
        !! The step limit.
    rst = this%m_stepLimit
end function

! --------------------
subroutine oi_set_step_limit(this, x)
    !! Sets the limit on the number of integration steps that may be taken 
    !! before the solver terminates.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    integer, intent(in) :: x
        !! The step limit.
    this%m_stepLimit = x
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_allow_overshoot(this) result(rst)
    !! Gets a value determining if the solver is allowed to overshoot the final
    !! value in the integration range.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    logical :: rst
        !! True if the solver can overshoot, and then interpolate to achieve the
        !! required final value; else, false thereby indicating the solver 
        !! cannot overshoot.
    rst = this%m_allowOvershoot
end function

! --------------------
subroutine oi_set_allow_overshoot(this, x)
    !! Sets a value determining if the solver is allowed to overshoot the final
    !! value in the integration range.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    logical, intent(in) :: x
        !! True if the solver can overshoot, and then interpolate to achieve the
        !! required final value; else, false thereby indicating the solver 
        !! cannot overshoot.
    this%m_allowOvershoot = x
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_log_failed_step(this, n)
    !! Logs a failed integration step.
    class(ode_integrator), intent(inout) :: this
        !! The [[ode_integrator]] object.
    integer, intent(in), optional :: n
        !! An optional input useful for logging more than one step.  If not
        !! supplied, this routine will log a single step.

    if (present(n)) then
        this%m_buffer%statistics%failed_steps = &
            this%m_buffer%statistics%failed_steps + n
    else
        this%m_buffer%statistics%failed_steps = &
            this%m_buffer%statistics%failed_steps + 1
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_log_successful_step(this, n)
    !! Logs a successful integration step.
    class(ode_integrator), intent(inout) :: this
        !! The [[ode_integrator]] object.
    integer, intent(in), optional :: n
        !! An optional input useful for logging more than one step.  If not
        !! supplied, this routine will log a single step.

    if (present(n)) then
        this%m_buffer%statistics%successful_steps = &
            this%m_buffer%statistics%successful_steps + n
    else
        this%m_buffer%statistics%successful_steps = &
            this%m_buffer%statistics%successful_steps + 1
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_log_jacobian_eval(this, n)
    !! Logs a Jacobian matrix evaluation.
    class(ode_integrator), intent(inout) :: this
        !! The [[ode_integrator]] object.
    integer, intent(in), optional :: n
        !! An optional input useful for logging more than one Jacobian matrix
        !! evaluation.  If not supplied, this routine will log a single
        !! Jacobian evaluation.

    if (present(n)) then
        this%m_buffer%statistics%jacobian_evaluations = &
            this%m_buffer%statistics%jacobian_evaluations + n
    else
        this%m_buffer%statistics%jacobian_evaluations = &
            this%m_buffer%statistics%jacobian_evaluations + 1
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_log_function_eval(this, n)
    !! Logs a function evaluation, or evaluations.
    class(ode_integrator), intent(inout) :: this
        !! The [[ode_integrator]] object.
    integer, intent(in), optional :: n
        !! An optional input useful for logging more than one function 
        !! evaluation.  If not supplied, this routine will log a single 
        !! function evaluation.
    if (present(n)) then
        this%m_buffer%statistics%function_evaluations = &
            this%m_buffer%statistics%function_evaluations + n
    else
        this%m_buffer%statistics%function_evaluations = &
            this%m_buffer%statistics%function_evaluations + 1
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_log_matrix_factorization(this, n)
    !! Logs a matrix factorization.
    class(ode_integrator), intent(inout) :: this
        !! The [[ode_integrator]] object.
    integer, intent(in), optional :: n
        !! An optional input useful for logging more than one matrix 
        !! factorization.  If not supplied, this routine will log a single
        !! factorization.
    if (present(n)) then
        this%m_buffer%statistics%matrix_factorizations = &
            this%m_buffer%statistics%matrix_factorizations + n
    else
        this%m_buffer%statistics%matrix_factorizations = &
            this%m_buffer%statistics%matrix_factorizations + 1
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_set_dense_enabled(this, enabled)
    !! Enables or disables explicit dense-output preparation.
    !! When enabled, the integrator prepares dense output for each accepted
    !! step independently of whether the requested output grid requires it.
    !! Disabling this setting does not prevent dense output needed to produce
    !! requested solution values. Any previously prepared step is invalidated.
    class(ode_integrator), intent(inout) :: this
        !! The ODE integrator.
    logical, intent(in) :: enabled
        !! True to enable explicit dense-output preparation; false to disable it.
    this%m_denseEnabled = enabled
    call this%invalidate_dense_output()
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_dense_enabled(this) result(rst)
    !! Gets the explicit dense-output preparation setting.
    class(ode_integrator), intent(in) :: this
        !! The ODE integrator.
    logical :: rst
        !! True if explicit dense-output preparation is enabled; otherwise,
        !! false. This setting is independent of preparation required to
        !! produce values on a requested output grid.
    rst = this%m_denseEnabled
end function

! ------------------------------------------------------------------------------
pure function oi_has_dense_output(this) result(rst)
    !! Determines whether a prepared accepted-step interpolant is available.
    class(ode_integrator), intent(in) :: this
        !! The ODE integrator.
    logical :: rst
        !! True if dense output is currently valid; otherwise, false.
    rst = this%m_denseValid
end function

! ------------------------------------------------------------------------------
subroutine oi_invalidate_dense(this)
    !! Marks any prepared dense-output interpolant as unavailable.
    class(ode_integrator), intent(inout) :: this
        !! The ODE integrator.
    this%m_denseValid = .false.
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_set_dense_interval(this, x_start, x_end, n)
    !! Publishes the interval covered by retained dense-output interpolation
    !! data for one accepted step. The endpoints are recorded in integration
    !! order, which may be increasing or decreasing. The interpolant is marked
    !! valid after the interval and state dimension have been recorded.
    !! An error stop occurs if either endpoint is non-finite, the endpoints
    !! are equal, or the state dimension is less than one.
    class(ode_integrator), intent(inout) :: this
        !! The ODE integrator.
    real(real64), intent(in) :: x_start
        !! The first endpoint of the accepted step in integration order.
    real(real64), intent(in) :: x_end
        !! The second endpoint of the accepted step in integration order.
    integer, intent(in) :: n
        !! The number of dependent variables represented by the interpolant.
    if (.not.ieee_is_finite(x_start) .or. .not.ieee_is_finite(x_end)) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    if (x_start == x_end .or. n < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_denseInterval = [x_start, x_end]
    this%m_denseSize = n
    this%m_denseValid = .true.
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_get_dense_interval(this, x_start, x_end)
    !! Gets the endpoints of the currently prepared dense-output interval.
    !! An error stop occurs if no valid interpolant is available.
    class(ode_integrator), intent(in) :: this
        !! The ODE integrator.
    real(real64), intent(out) :: x_start
        !! Receives the first endpoint in integration order.
    real(real64), intent(out) :: x_end
        !! Receives the second endpoint in integration order.
    if (.not.this%m_denseValid) error stop DIFFEQ_INVALID_OPERATION_ERROR
    x_start = this%m_denseInterval(1)
    x_end = this%m_denseInterval(2)
end subroutine

! ------------------------------------------------------------------------------
subroutine oi_evaluate_dense(this, x, y)
    !! Evaluates the interpolant for the last prepared accepted step at a
    !! requested independent-variable value. The query must be finite and
    !! within the step interval, including either endpoint. This operation
    !! does not advance or otherwise change solver state and does not evaluate
    !! the ODE function. An error stop occurs if no valid interpolant is
    !! available, the output array has the wrong size, or the query is
    !! non-finite or outside the interval.
    class(ode_integrator), intent(in) :: this
        !! The ODE integrator containing the prepared interpolant.
    real(real64), intent(in) :: x
        !! The independent-variable value at which to evaluate the solution.
    real(real64), intent(out) :: y(:)
        !! Receives the interpolated dependent-variable values. Its size must
        !! equal the number of dependent variables in the prepared step.
    if (.not.this%m_denseValid) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (size(y) /= this%m_denseSize) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (.not.ieee_is_finite(x)) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (x < minval(this%m_denseInterval) .or. x > maxval(this%m_denseInterval)) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    call this%compute_dense_output(x, y)
end subroutine

! ******************************************************************************
! ODE_EVENT_SOLVER_ROUTINES
! ------------------------------------------------------------------------------
pure function oes_get_max_iter(this) result(rst)
    !! Gets the maximum number of iterations allowed for the root-finding solver.
    class(ode_event_solver), intent(in) :: this
        !! The [[ode_event_solver]] object.
    integer :: rst
        !! The iteration count limit.
    rst = this%m_maxIter
end function

! --------------------
subroutine oes_set_max_iter(this, x)
    !! Sets the maximum number of iterations allowed for the root-finding solver.
    class(ode_event_solver), intent(inout) :: this
        !! The [[ode_event_solver]] object.
    integer, intent(in) :: x
        !! The iteration count limit.
    this%m_maxIter = x
end subroutine

! ------------------------------------------------------------------------------
pure function oes_get_tol(this) result(rst)
    !! Gets the root-finding solver tolerance.
    class(ode_event_solver), intent(in) :: this
        !! The [[ode_event_solver]] object.
    real(real64) :: rst
        !! The tolerance value.
    rst = this%m_tol
end function

! --------------------
subroutine oes_set_tol(this, x)
    !! Sets the root-finding solver tolerance.
    class(ode_event_solver), intent(inout) :: this
        !! The [[ode_event_solver]] object.
    real(real64), intent(in) :: x
        !! The tolerance value.
    this%m_tol = x
end subroutine

! ------------------------------------------------------------------------------
function oes_find_event(this, fcn, integrator, neqn, x1, y1, x2, y2, x, y, &
    args) result(rst)
    !! Attempts to find an event via Brent's method for root finding for each
    !! of the M event equations (M does not have to be equal to the number
    !! of state equations).  If no events are found, this function simply
    !! returns false.  If a root, or roots, are found, the first occurring
    !! root is returned, and iteration statistics are updated in the
    !! [[ode_event_solver]].
    !!
    !! If no event function has been defined a DIFFEQ_NULL_POINTER_ERROR is
    !! thrown.
    class(ode_event_solver), intent(inout) :: this
        !! The [[ode_event_solver]] object.
    procedure(ode_event_function), intent(in), pointer :: fcn
        !! The [[ode_event_function]].
    class(ode_integrator), intent(in) :: integrator
        !! The [[ode_integrator]] being used to integrate the ODE's.
    integer, intent(in) :: neqn
        !! The number of event equations.  This value must be at least 1; else,
        !! a DIFFEQ_INVALID_INPUT_ERROR is thrown.
    real(real64), intent(in) :: x1
        !! The initial bound of the search region.
    real(real64), intent(in), dimension(:) :: y1
        !! The N-element state vector at x1.
    real(real64), intent(in) :: x2
        !! The terminal bound of the search region.  This value must be
        !! different than x2 by an amount larger than machine epsilon; else,
        !! a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    real(real64), intent(in), dimension(:) :: y2
        !! The N-element state vector at x2.  This array must be sized to match
        !! x1; else, a DIFFEQ_ARRAY_SIZE_ERROR is thrown.
    real(real64), intent(out) :: x
        !! The value of the independent variable at which the root is found.
        !! If no root is found, this value is returned as a NaN.
    real(real64), intent(out), dimension(:) :: y
        !! The N-element array containing the state vector at x.  If no root
        !! is found, this array is populated with NaN.  This array must be
        !! sized to match x1; else, a DIFFEQ_ARRAY_SIZE_ERROR is thrown.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information in and out
        !! of the routine.
    logical :: rst
        !! True if the process converged to a root; else, false.

    ! Local Variables
    integer :: i, nstates, niter, neval, ind
    real(real64) :: eps
    real(real64), allocatable, dimension(:) :: xs
    real(real64), allocatable, dimension(:,:) :: ys
    logical, allocatable, dimension(:) :: found_event

    ! Input Validation
    nstates = size(y1)
    eps = epsilon(eps)
    if (size(y2) /= nstates .or. size(y) /= nstates) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (abs(x2 - x1) <= eps) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (neqn < 1) error stop DIFFEQ_INVALID_INPUT_ERROR

    ! Initialization
    allocate( &
        xs(neqn), &
        ys(nstates, neqn), &
        found_event(neqn) &
    )
    x = ieee_value(0.0d0, ieee_quiet_nan)
    y = ieee_value(0.0d0, ieee_quiet_nan)
    rst = .false.

    ! Reset the iteration statistics information
    this%solver_statistics%event_index = 0
    this%solver_statistics%function_evaluations = 0
    this%solver_statistics%interation_count = 0

    ! Locate events
    do i = 1, neqn
        found_event(i) = this%find_channel_event(fcn, integrator, i, neqn, &
            x1, y1, x2, y2, xs(i), ys(:,i), niter, neval, args = args)
    end do
    rst = any(found_event)
    if (.not.rst) return

    ! If any events were found, find the first event
    ind = 0
    do i = 1, neqn
        if (.not.found_event(i)) cycle
        if (ind == 0) then
            ind = i
        else if (x1 < x2 .and. xs(i) < xs(ind)) then
            ind = i
        else if (x1 > x2 .and. xs(i) > xs(ind)) then
            ind = i
        end if
    end do
    x = xs(ind)
    y = ys(:,ind)

    ! Store the iteration statistics information
    this%solver_statistics%function_evaluations = neval
    this%solver_statistics%interation_count = niter
    this%solver_statistics%event_index = ind
end function

! ------------------------------------------------------------------------------
function oes_find_single_channel_event(this, fcn, integrator, index, nevents, &
    x1, y1, x2, y2, x, y, niter, neval, args) result(rst)
    !! Performs the root finding iteration process via Brent's method for a
    !! single channel.  As this routine is intended to stay private, there is
    !! no error checking that is performed.  However, there is a call to other
    !! type-bound procedures that do check their inputs courtesy of their public
    !! interface; therefore, it is still possible that this routine would 
    !! trigger an 'error stop' condition.  For this reason, be sure of the 
    !! inputs supplied to this routine.
    class(ode_event_solver), intent(in) :: this
        !! The [[ode_event_solver]] object.
    procedure(ode_event_function), intent(in), pointer :: fcn
        !! The [[ode_event_function]].
    class(ode_integrator), intent(in) :: integrator
        !! The [[ode_integrator]] being used to integrate the ODE's.
    integer, intent(in) :: index
        !! The index of the channel on which to operate.
    integer, intent(in) :: nevents
        !! The number of event equations.
    real(real64), intent(in) :: x1
        !! The initial bound of the search region.
    real(real64), intent(in), dimension(:) :: y1
        !! The N-element state vector at x1.
    real(real64), intent(in) :: x2
        !! The terminal bound of the search region.
    real(real64), intent(in), dimension(:) :: y2
        !! The N-element state vector at x2.
    real(real64), intent(out) :: x
        !! The value of the independent variable at which the root is found.
        !! If no root is found, this value is returned as a NaN.
    real(real64), intent(out), dimension(:) :: y
        !! The N-element array containing the state vector at x.  If no root
        !! is found, this array is populated with NaN.
    integer, intent(out) :: niter
        !! The number of iterations taken.
    integer, intent(out) :: neval
        !! The number of event function evaluations taken.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information in and out
        !! of the routine.
    logical :: rst
        !! True if the process converged to a root; else, false.

    ! Local Variables
    integer :: maxiter, nstates
    real(real64), allocatable, dimension(:) :: g1, g2, g, yinterp
    real(real64) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, eps, tol, tol1, xm, &
        min1, min2, mn

    ! Initialization
    rst = .false.
    x = ieee_value(0.0d0, ieee_quiet_nan)
    y = ieee_value(0.0d0, ieee_quiet_nan)
    niter = 0
    maxiter = this%get_max_allowed_iterations()
    eps = epsilon(eps)
    tol = this%get_solver_tolerance()
    nstates = size(y)
    allocate( &
        g1(nevents), &
        g2(nevents), &
        g(nevents), &
        yinterp(nstates) &
    )

    ! Evaluate the event functions.
    call fcn(x1, y1, g1, args = args)
    call fcn(x2, y2, g2, args = args)
    neval = 2

    ! Test to see if a root exists within this range for this index
    if (.not.this%is_event_suspected(g1, g2, index)) then
        ! No root is expected in this bracket, for this channel - return
        rst = .false.
        return
    end if

    ! Iteration Process
    a = x1
    b = x2
    c = x2
    fa = g1(index)
    fb = g2(index)
    fc = fb
    do niter = 1, maxiter
        ! Adjust the bounding interval
        if ((fb > 0.0d0 .and. fc > 0.0d0) .or. &
            (fb < 0.0d0 .and. fc < 0.0d0)) &
        then
            c = a
            fc = fa
            e = b - a
            d = e
        end if
        if (abs(fc) < abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        end if

        ! Convergence Check
        tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
        xm = 0.5d0 * (c - b)
        if (abs(xm) <= tol1 .or. fb == 0.0d0) then
            ! b is the root
            rst = .true.
            x = b
            call integrator%evaluate_dense(x, y)
            exit
        end if

        if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
            ! Attempt interpolation
            s = fb / fa
            if (a == c) then
                p = 2.0d0 * xm * s
                q = 1.0d0 - s
            else
                q = fa / fc
                r = fb / fc
                p = s * (2.0d0 * xm * q * (q - r) - (b - a) * (r - 1.0d0))
                q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
            end if

            ! Bounds checking
            if (p > 0.0d0) q = -q
            p = abs(p)
            min1 = 3.0d0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)
            if (min1 < min2) then
                mn = min1
            else
                mn = min2
            end if
            if (2.0d0 * p < mn) then
                ! Accept the interpolation
                e = d
                d = p / q
            else
                ! The interpolation failed, use bisection
                d = xm
                e = d
            end if
        else
            ! The bounds are decreasing too slowly, use bisection
            d = xm
            e = d
        end if

        ! Update
        a = b
        fa = fb

        ! Evaluate the new trial root
        if (abs(d) > tol1) then
            b = b + d
        else
            b = b + sign(tol1, xm)
        end if
        call integrator%evaluate_dense(b, yinterp)
        call fcn(b, yinterp, g, args = args)
        fb = g(index)
        neval = neval + 1
    end do
end function

! ------------------------------------------------------------------------------
pure function oes_event_in_bracket(this, f1, f2, index) result(rst)
    !! Tests to see if a root is likely to exist in the given bracket.
    class(ode_event_solver), intent(in) :: this
        !! The [[ode_event_solver]] object.
    real(real64), intent(in), dimension(:) :: f1
        !! An M-element array containing the values of the M event functions at
        !! the initial interval point.
    real(real64), intent(in), dimension(:) :: f2
        !! An M-element array containing the values of the M event functions at
        !! the terminal interval point.  An DIFFEQ_ARRAY_SIZE_ERROR is thrown
        !! if this array is not the same size as f1.
    integer(int32), intent(in) :: index
        !! The index of the event function to validate.  The index must exist
        !! on the interval [1, M], where M is the number of event functions.  A
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown if index is outside this 
        !! interval.
    logical :: rst
        !! True if an event is suspected to exist within the bracket; else,
        !! false.

    ! Input Checking
    if (size(f1) /= size(f2)) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (index < 1 .or. index > size(f1)) error stop DIFFEQ_INVALID_INPUT_ERROR

    ! Process
    rst = f1(index) * f2(index) <= 0.0d0
end function

! ------------------------------------------------------------------------------
end module