module diffeq_base
    !! A collection of base types for the DIFFEQ library.
    use iso_fortran_env
    use diffeq_errors
    use ferror
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
    public :: single_step_pre_step_routine
    public :: single_step_interpolate
    public :: single_step_integer_inquiry
    public :: single_step_integrator

! ------------------------------------------------------------------------------
    interface
        subroutine ode(x, y, dydx)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        subroutine ode_jacobian(x, y, jac)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        subroutine ode_mass_matrix(x, y, m)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: m
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type ode_container
        !! A container for the routine containing the ODEs to integrate.
        logical, private :: m_massDependent = .true.
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
    end type

! ------------------------------------------------------------------------------
    type, abstract :: ode_integrator
        !! The most basic ODE integrator object capable of integrating
        !! systems of ODE's.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            ! The internal solution buffer.
        integer(int32), private :: m_bufferCount = 0
            ! The number of solution points stored in the buffer.
        real(real64), private :: m_abstol = 1.0d-6
            ! The absolute tolerance value applied to each equation.
        real(real64), private :: m_reltol = 1.0d-6
            ! The relative tolerance value applied to each equation.
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
        integer(int32), private :: m_stepLimit = 1000000
            ! A limit on the total number of integration steps.  Exceeding this
            ! value may indicate that a different integrator should be
            ! considered.  Either that, or this is a really large-sized
            ! problem that is being solved.
        logical, private :: m_allowOvershoot = .true.
            ! True if the solver is allowed to overshoot the final value and
            ! interpolate back.  False, if the solver must terminate on the
            ! final value.
    contains
        procedure(ode_solver), public, pass, deferred :: solve
            !! Solves the supplied system of ODE's.
        procedure(ode_integer_inquiry), public, pass, deferred :: get_order
            !! Returns the order of the integrator.
        procedure, public :: append_to_buffer => oi_append_to_buffer
            !! Appends the supplied solution point to the internal solution 
            !! buffer.
        procedure, public :: get_solution => oi_get_solution
            !! Returns the solution computed by the integrator.
        procedure, public :: clear_buffer => oi_clear_buffer
            !! Clears the contents of the buffer.
        procedure, public :: get_absolute_tolerance => oi_get_abs_tol
            !! Gets the absolute error tolerance.
        procedure, public :: set_absolute_tolerance => oi_set_abs_tol
            !! Sets the absolute error tolerance.
        procedure, public :: get_relative_tolerance => oi_get_abs_tol
            !! Gets the relative error tolerance.
        procedure, public :: set_relative_tolerance => oi_set_abs_tol
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
    end type

    interface
        subroutine ode_solver(this, sys, x, iv, err)
            !! Solves the supplied system of ODE's.
            use iso_fortran_env
            use ferror
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
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to 
                !! provide error handling.  Possible errors and warning messages
                !! that may be encountered are as follows.
                !!
                !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
                !!      memory allocation issue.
                !!
                !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is 
                !!      defined.
                !!
                !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if there are less than 
                !!      2 values given in the independent variable array x.
        end subroutine

        pure function ode_integer_inquiry(this) result(rst)
            !! Returns an integer value from the ode_integrator object.
            use iso_fortran_env
            import ode_integrator
            class(ode_integrator), intent(in) :: this
                !! The ode_integrator object.
            integer(int32) :: rst
                !! The requested value.
        end function
    end interface

! ------------------------------------------------------------------------------
    type, abstract, extends(ode_integrator) :: single_step_integrator
        !! The most basic, single-step integrator object capable of integrating
        !! systems of ODE's.
    contains
        procedure(attempt_single_step), public, pass, deferred :: attempt_step
            !! Attempts an integration step for a single-step integrator.
        procedure(get_single_step_logical_parameter), public, pass, &
            deferred :: get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure(single_step_post_step_routine), public, pass, deferred :: &
            post_step_action
            !! Performs actions such as setting up interpolation after 
            !! completion of a successful integration step.
        procedure(single_step_interpolate), public, pass, deferred :: &
            interpolate
            !! Performs an interpolation to estimate the solution at the
            !! requested point.
        procedure(single_step_pre_step_routine), public, pass, deferred :: &
            pre_step_action
            !! Provides a routine for performing any actions, such as setting
            !! up Jacobian calculations.
        procedure(single_step_integer_inquiry), public, pass, deferred :: &
            get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure, public :: solve => ssi_ode_solver
            !! Solves the supplied system of ODE's.
    end type

    interface
        subroutine attempt_single_step(this, sys, h, x, y, f, yn, fn, yerr, k)
            use iso_fortran_env
            import single_step_integrator
            import ode_container
            !! Attempts an integration step for a single-step integrator.
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current solution at x.
            real(real64), intent(in), dimension(:) :: f
                !! An N-element array containing the values of the derivatives
                !! at x.
            real(real64), intent(out), dimension(:) :: yn
                !! An N-element array where this routine will write the next
                !! solution estimate at x + h.
            real(real64), intent(out), dimension(:) :: fn
                !! An N-element array where this routine will write the next
                !! derivative estimate at x + h.
            real(real64), intent(out), dimension(:) :: yerr
                !! An N-element array where this routine will write an estimate
                !! of the error in each equation.
            real(real64), intent(out), dimension(:,:) :: k
                !! An N-by-NSTAGES matrix containing the derivatives at each
                !! stage.
        end subroutine

        pure function get_single_step_logical_parameter(this) result(rst)
            !! Returns a logical parameter from a single_step_integrator object.
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
            logical :: rst
                !! The parameter.
        end function

        subroutine single_step_post_step_routine(this, sys, dense, x, xn, y, &
            yn, f, fn, k)
            !! Provides a routine for performing any actions, such as setting
            !! up interpolation, after successful completion of a step.
            use iso_fortran_env
            import single_step_integrator
            import ode_container
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            logical, intent(in) :: dense
                !! Determines if dense output is requested (true); else, false.
            real(real64), intent(in) :: x
                !! The previous value of the independent variable.
            real(real64), intent(in) :: xn
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the solution at x.
            real(real64), intent(in), dimension(:) :: yn
                !! An N-element array containing the solution at xn.
            real(real64), intent(in), dimension(:) :: f
                !! An N-element array containing the derivatives at x.
            real(real64), intent(in), dimension(:) :: fn
                !! An N-element array containing the derivatives at xn.
            real(real64), intent(inout), dimension(:,:) :: k
                !! An N-by-NSTAGES matrix containing the derivatives at each
                !! stage.
        end subroutine

        subroutine single_step_interpolate(this, x, xn, yn, fn, xn1, yn1, &
            fn1, y)
            !! Provides a routine for interpolation.
            use iso_fortran_env
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
            real(real64), intent(in) :: x
                !! The value of the independent variable at which to compute
                !! the interpolation.
            real(real64), intent(in) :: xn
                !! The previous value of the independent variable at which the
                !! solution is computed.
            real(real64), intent(in), dimension(:) :: yn
                !! An N-element array containing the solution at xn.
            real(real64), intent(in), dimension(:) :: fn
                !! An N-element array containing the derivatives at xn.
            real(real64), intent(in) :: xn1
                !! The value of the independent variable at xn + h.
            real(real64), intent(in), dimension(:) :: yn1
                !! An N-element array containing the solution at xn + h.
            real(real64), intent(in), dimension(:) :: fn1
                !! An N-element array containing the derivatives at xn + h.
            real(real64), intent(out), dimension(:) :: y
                !! An N-element array where this routine will write the 
                !! solution values interpolated at x.
        end subroutine

        subroutine single_step_pre_step_routine(this, prevs, sys, h, x, y, f, &
            err)
            !! Provides a routine for performing any actions, such as setting
            !! up Jacobian calculations.
            use iso_fortran_env
            use ferror
            import single_step_integrator
            import ode_container
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
            logical, intent(in) :: prevs
                !! Defines the status of the previous step.  The value is true 
                !! if the previous step was successful; else, false if the 
                !! previous step failed.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current solution at x.
            real(real64), intent(in), dimension(:) :: f
                !! An N-element array containing the values of the derivatives
                !! at x.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to 
                !! provide error handling.
        end subroutine

        pure function single_step_integer_inquiry(this) result(rst)
            !! Gets an integer from the integrator.
            use iso_fortran_env
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
            integer(int32) :: rst
                !! The integer value.
        end function
    end interface
  
! ------------------------------------------------------------------------------    ! TO DO: multi-step integrator

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
subroutine oc_jacobian(this, x, y, jac, err)
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling. Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined,
        !!      and the calculation is being performed by finite differences.
        !! 
        !!  - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if jac is not N-by-N.

    ! Local Variables
    integer(int32) :: i, ndof
    real(real64) :: h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    ndof = size(y)
    h = this%get_finite_difference_step()

    ! Input Checking
    if (size(jac, 1) /= ndof .or. size(jac, 2) /= ndof) then
        call report_matrix_size_error(errmgr, "oc_jacobian", "jac", &
            ndof, ndof, size(jac, 1), size(jac, 2))
        return
    end if

    ! Use a user-defined routine, and then be done
    if (associated(this%jacobian)) then
        call this%jacobian(x, y, jac)
        return
    end if

    ! Allocate workspace.  No action is taken if the proper workspace is
    ! already allocated.
    call this%allocate_workspace(ndof, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Finite Difference Approximation
    ! J(i,j) = df(i) / dy(j)
    this%m_jwork(1:ndof) = y
    call this%fcn(x, y, this%m_jwork(ndof+1:))
    do i = 1, ndof
        this%m_jwork(i) = this%m_jwork(i) + h
        call this%fcn(x, this%m_jwork(1:ndof), jac(:,i))
        jac(:,i) = (jac(:,i) - this%m_jwork(ndof+1:)) / h
        this%m_jwork(i) = y(i)
    end do

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine oc_alloc_workspace(this, ndof, err)
    ! Use to allocate internal workspaces.  This routine only takes action
    ! if the workspace array(s) are not sized properly for the application.
    class(ode_container), intent(inout) :: this
        ! The ode_container object.
    integer(int32), intent(in) :: ndof
        ! The number of degrees of freedom.
    class(errors), intent(inout) :: err
        ! The error handling object.

    ! Local Variables
    integer(int32) :: flag

    ! Jacobian Workspace Allocation
    if (allocated(this%m_jwork)) then
        if (size(this%m_jwork) /= 2 * ndof) then
            deallocate(this%m_jwork)
            allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
            if (flag /= 0) then
                call report_memory_error(err, "oc_alloc_workspace", flag)
                return
            end if
        end if
    else
        allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
        if (flag /= 0) then
            call report_memory_error(err, "oc_alloc_workspace", flag)
            return
        end if
    end if

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
pure function oc_get_is_ode_defined(this) result(rst)
    !! Gets a logical value determining if the ODE routine has been defined.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    logical :: rst
        !! True if the ODE routine has been defined; else, false.

    rst = associated(this%fcn)
end function

! ******************************************************************************
! ODE_INTEGRATOR ROUTINES
! ------------------------------------------------------------------------------
subroutine oi_append_to_buffer(this, x, y, err)
    !! Appends the supplied solution point to the internal solution buffer.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The independent variable value.
    real(real64), intent(in), dimension(:) :: y
        !! The values of the dependent variables corresponding to x.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
        !!      memory allocation issue.

    ! Parameters
    integer(int32), parameter :: buffer_size = 1000

    ! Local Variables
    integer(int32) :: i, start, m, n, neqn, flag
    real(real64), allocatable, dimension(:,:) :: copy
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (allocated(this%m_buffer)) then
        neqn = size(this%m_buffer, 2) - 1
    else
        neqn = size(y)
    end if
    n = neqn + 1
    start = this%m_bufferCount + 1

    ! Input Checking
    if (size(y) /= neqn) then
        call report_array_size_error(errmgr, "oi_append_to_buffer", "y", neqn, &
            size(y))
        return
    end if

    ! Allocate memory if necessary
    if (.not.allocated(this%m_buffer)) then
        allocate(this%m_buffer(buffer_size, n), stat = flag)
        if (flag /= 0) go to 10
    else
        m = size(this%m_buffer, 1)
        if (start == m + 1) then
            allocate(copy(m, n), source = this%m_buffer, stat = flag)
            if (flag /= 0) go to 10
            m = m + buffer_size
            deallocate(this%m_buffer)
            allocate(this%m_buffer(m, n), stat = flag)
            if (flag /= 0) go to 10
            this%m_buffer(:this%m_bufferCount,:) = copy
        end if
    end if

    ! Store the result
    this%m_buffer(start,1) = x
    this%m_buffer(start,2:) = y
    this%m_bufferCount = this%m_bufferCount + 1

    ! End
    return

    ! Memory Error Handler
10  continue
    call report_memory_error(errmgr, "oi_append_to_buffer", flag)
    return
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

    if (allocated(this%m_buffer)) then
        rst = this%m_buffer(:this%m_bufferCount,:)
    else
        allocate(rst(0, 0))
    end if
end function

! ------------------------------------------------------------------------------
subroutine oi_clear_buffer(this)
    !! Clears the contents of the buffer.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.

    ! Clear the buffer
    if (allocated(this%m_buffer)) deallocate(this%m_buffer)
    this%m_bufferCount = 0
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_abs_tol(this) result(rst)
    !! Gets the absolute error tolerance.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The tolerance value.
    rst = this%m_abstol
end function

! --------------------
subroutine oi_set_abs_tol(this, x)
    !! Sets the absolute error tolerance.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance value.
    this%m_abstol = x
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_rel_tol(this) result(rst)
    !! Gets the relative error tolerance.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    real(real64) :: rst
        !! The tolerance value.
    rst = this%m_reltol
end function

! --------------------
subroutine oi_set_rel_tol(this, x)
    !! Sets the relative error tolerance.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance value.
    this%m_reltol = x
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
    integer(int32) :: i, n
    real(real64) :: sf, atol, rtol

    ! Initialization
    n = size(y)
    atol = this%get_absolute_tolerance()
    rtol = this%get_relative_tolerance()

    ! Process
    rst = 0.0d0
    do i = 1, n
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
function oi_next_step(this, e, eold, h, x, err) result(rst)
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.
        !!
        !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size
        !!      becomes too small in magnitude.
    real(real64) :: rst
        !! The new step size estimate.

    ! Parameters
    real(real64), parameter :: minscale = 2.0d-1
    real(real64), parameter :: maxscale = 1.0d1

    ! Local Variables
    integer(int32) :: k
    real(real64) :: alpha, beta, fs, maxstep, minstep, scale
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
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
    if (abs(rst) < minstep) then
        call report_step_size_too_small(errmgr, "oi_next_step", x, rst)
        return
    end if
end function

! ------------------------------------------------------------------------------
subroutine oi_initial_step(this, sys, xo, xf, yo, fo, h)
    !! Computes an estimate of an initial step size.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    class(ode_container), intent(in) :: sys
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

    ! Local Variables
    real(real64) :: e, dx

    ! Use a very basic estimate of an initial step size.  The catch is that a
    ! single function evaluation must be made; however, this is likely needed
    ! in the first place, so no real extra work is necessary.
    dx = 0.1d0 * (xf - xo)
    e = max(this%get_absolute_tolerance(), this%get_relative_tolerance())
    call sys%fcn(xo, yo, fo)
    h = 2.0d0 * e / norm2(fo)
    h = min(abs(dx), h)
    h = sign(h, dx)
end subroutine

! ------------------------------------------------------------------------------
pure function oi_get_step_limit(this) result(rst)
    !! Gets the limit on the number of integration steps that may be taken 
    !! before the solver terminates.
    class(ode_integrator), intent(in) :: this
        !! The ode_integrator object.
    integer(int32) :: rst
        !! The step limit.
    rst = this%m_stepLimit
end function

! --------------------
subroutine oi_set_step_limit(this, x)
    !! Sets the limit on the number of integration steps that may be taken 
    !! before the solver terminates.
    class(ode_integrator), intent(inout) :: this
        !! The ode_integrator object.
    integer(int32), intent(in) :: x
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

! ******************************************************************************
! SINGLE-STEP INTEGRATOR
! ------------------------------------------------------------------------------
subroutine ssi_ode_solver(this, sys, x, iv, err)
    !! Solves the supplied system of ODE's.
    class(single_step_integrator), intent(inout) :: this
        !! The single_step_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array of independent variable values at which to return the 
        !! the solution to the ODE's.
    real(real64), intent(in), dimension(:) :: iv
        !! An array containing the initial values for each ODE.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
        !!      memory allocation issue.
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is 
        !!      defined.
        !!
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if there are less than 
        !!      2 values given in the independent variable array x.

    ! Local Variables
    logical :: dense, success
    integer(int32) :: i, j, n, neqn, flag, nsteps, nstages
    real(real64) :: h, xo, xn, xmax, ei, eold
    real(real64), allocatable, dimension(:) :: f, y, yn, fn, yerr, yi
    real(real64), allocatable, dimension(:,:) :: k
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    success = .true.

    ! Input Checking
    if (n < 2) then
        call report_array_size_error(errmgr, "ssi_ode_solver", "x", 2, n)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "ssi_ode_solver")
        return
    end if

    ! Additional Initialization
    neqn = size(iv)
    xo = x(1)
    xmax = x(n)
    j = 2
    eold = 1.0d-4
    dense = (n > 2)
    nsteps = this%get_step_limit()
    nstages = this%get_stage_count()
    
    ! Memory Allocations
    allocate( &
        f(neqn), &
        y(neqn), &
        yn(neqn), &
        fn(neqn), &
        yerr(neqn), &
        k(neqn, nstages),  &
        stat = flag &
    )
    if (flag == 0 .and. dense) allocate(yi(neqn), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "ssi_ode_solver", flag)
        return
    end if

    ! Estimate an initial step size
    !
    ! Outputs:
    ! - f: Value of the derivatives at xo
    ! - h: Initial step size estimate
    call this%estimate_inital_step_size(sys, xo, xmax, iv, f, h)
    
    ! Store the initial conditions
    call this%append_to_buffer(x(1), iv, errmgr)
    if (errmgr%has_error_occurred()) return
    y = iv

    ! Cycle until integration is complete
    do i = 1, nsteps
        ! Perform any pre-step actions
        call this%pre_step_action(success, sys, h, xo, y, f, errmgr)
        if (errmgr%has_error_occurred()) return
        
        ! Attempt a step
        call this%attempt_step(sys, h, xo, y, f, yn, fn, yerr, k)
        xn = xo + h

        ! Compute a normalized error value.  A value < 1 indicates success
        ei = this%compute_error_norm(y, yn, yerr)

        ! Determine the next step size
        h = this%estimate_next_step_size(ei, eold, h, xo, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Reject the step?
        success = ei <= 1.0d0
        if (.not.success) cycle   ! We failed.  Try again with a new step size

        ! If we're here, the step has been successful.  Take any post-step
        ! action such as setting up interpolation routines, etc.
        call this%post_step_action(sys, dense, xo, xn, y, yn, f, fn, k)

        ! Do we need to interpolate for dense output, or can we just store
        ! values and move on
        if (dense) then
            ! Perform the interpolation as needed until a new step is required
            interp : do while (abs(x(j)) <= abs(xn))
                call this%interpolate(x(j), xo, y, f, xn, yn, fn, yi)
                call this%append_to_buffer(x(j), yi, errmgr)
                if (errmgr%has_error_occurred()) return
                j = j + 1
                if (j > n) exit interp
            end do interp
        else
            ! Store the values and move on
            call this%append_to_buffer(xn, yn, errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Are we done?
        if (abs(xn) >= abs(xmax)) then
            ! Deal with the case where output is only returned at the
            ! integration points and the solver oversteps the endpoint.
            if (abs(xn) > abs(xmax)) then
                ! Interpolate to get the solution at xmax
                call this%post_step_action(sys, .true., xo, xn, y, yn, f, fn, k)
                call this%interpolate(xmax, xo, y, f, xn, yn, fn, yi)
                call this%append_to_buffer(xmax, yi, errmgr)
                if (errmgr%has_error_occurred()) return
            end if
            
            ! We're done
            go to 100
        end if

        ! Do we need to limit the step size to not overshoot the terminal value?
        if (this%get_allow_overshoot() .and. abs(xn + h) > abs(xmax)) then
            ! Limit the step size
            h = xmax - xn
        end if

        ! Update parameters
        xo = xn
        y = yn

        ! Is this an FSAL integrator?  If so, we already have the derivative
        ! values.  If not, we need to recompute the derivative values for the
        ! next iteration
        if (this%get_is_fsal()) then
            ! We have the derivative estimate already
            f = fn
        else
            ! Update the derivative estimates
            call sys%fcn(xn, yn, f) ! write to f, not fn
        end if
    end do

    ! If we're here, the solver has run out of allowable steps
    call report_excessive_integration_steps(errmgr, "ssi_ode_solver", nsteps, &
        xn)
    return

    ! End
100 continue
    return
end subroutine

! ------------------------------------------------------------------------------
end module