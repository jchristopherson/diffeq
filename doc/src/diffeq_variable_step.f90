module diffeq_variable_step
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    implicit none
    private
    public :: variable_step_integrator
    public :: variable_step_attempt
    public :: variable_step_action
    public :: variable_step_interpolation
    public :: next_step_size_calculator

    type, abstract, extends(ode_integrator) :: variable_step_integrator
        !! Defines a variable-step integrator.
        real(real64), private :: m_safetyfactor = 0.9d0
        real(real64), private :: m_maxstep = huge(1.0d0)
        real(real64), private :: m_minstep = 1.0d2 * epsilon(1.0d0)
        integer(int32), private :: m_maxitercount = 1000
        integer(int32), private :: m_maxstepcount = 1000000
        real(real64), private :: m_stepSize = 1.0d0
        real(real64), private :: m_nextStep = 1.0d0
        real(real64), private :: m_enormPrev = 1.0d0
        logical, private :: m_respectXMax = .true.
        ! Solution buffer
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
        integer(int32) :: m_bufferCount = 0
        ! Workspaces
        real(real64), private, allocatable, dimension(:) :: m_ework ! NEQN element
        ! Tolerances
        real(real64), private, allocatable, dimension(:) :: m_rtol ! NEQN element
        real(real64), private, allocatable, dimension(:) :: m_atol ! NEQN element
    contains
        procedure, private :: initialize => vsi_alloc_workspace
            !! Initializes the integrator.
        procedure(variable_step_attempt), deferred, public :: attempt_step
            !! Attempts a single integration step.
        procedure(variable_step_action), deferred, public :: on_successful_step
            !! Perform necessary actions on completion of a successful step.
        procedure, public :: get_safety_factor => vsi_get_safety_factor
            !! Gets a safety factor used to limit the predicted step size.
        procedure, public :: set_safety_factor => vsi_set_safety_factor
            !! Sets a safety factor used to limit the predicted step size.
        procedure, public :: get_max_step_size => vsi_get_max_step
            !! Gets the maximum allowed step size.
        procedure, public :: set_max_step_size => vsi_set_max_step
            !! Sets the maximum allowed step size.
        procedure, public :: get_min_step_size => vsi_get_min_step
            !! Gets the minimum allowed step size.
        procedure, public :: set_min_step_size => vsi_set_min_step
            !! Sets the minimum allowed step size.
        procedure, public :: get_max_per_step_iteration_count => &
            vsi_get_max_iter_count
            !! Gets the maximum number of iterations per step allowed.
        procedure, public :: set_max_per_step_iteration_count => &
            vsi_set_max_iter_count
            !! Sets the maximum number of iterations per step allowed.
        procedure, public :: get_max_integration_step_count => & 
            vsi_get_max_step_count
            !! Gets the maximum number of integration steps allowed.
        procedure, public :: set_max_integration_step_count => &
            vsi_set_max_step_count
            !! Sets the maximum number of integration steps allowed.
        procedure(next_step_size_calculator), deferred, public :: &
            compute_next_step_size
            !! Computes the next step size.
        procedure, public :: buffer_results => vsi_append_to_buffer
            !! Buffers a results set.
        procedure, public :: get_buffer_size => vsi_get_buffer_count
            !! Gets the number of entries into the solution buffer.
        procedure, public :: clear_buffer => vsi_clear_buffer
            !! Clears the results buffer.
        procedure, public :: get_step_size => vsi_get_step_size
            !! Gets the current step size.
        procedure, public :: set_step_size => vsi_set_step_size
            !! Sets the current step size.
        procedure, public :: get_next_step_size => vsi_get_next_step_size
            !! Gets the next step size.
        procedure, public :: set_next_step_size => vsi_set_next_step_size
            !! Sets the next step size.
        procedure, public :: get_respect_x_max => vsi_get_respect_xmax
            !! Gets a value determining if the integrator should respect a
            !! hard limit in the independent variable range.  If false, the 
            !! integrator may step pass the limit.
        procedure, public :: set_respect_x_max => vsi_set_respect_xmax
            !! Sets a value determining if the integrator should respect a
            !! hard limit in the independent variable range.  If false, the 
            !! integrator may step pass the limit.
        procedure, public :: step => vsi_step
            !! Takes one integration step.
        procedure, public :: estimate_first_step_size => vsi_estimate_first_step
            !! Computes an estimate to the first step size based upon the 
            !! initial function values.
        procedure(variable_step_interpolation), public, deferred :: interpolate
            !! Provides interpolation between integration points.
        procedure, public :: get_default_relative_tolerance => &
            vsi_get_default_rel_tol
            !! Gets the default relative error tolerance.
        procedure, public :: get_default_absolute_tolerance => &
            vsi_get_default_abs_tol
            !! Gets the default absolute error tolerance.
        procedure, public :: get_buffer_contents => vsi_get_buffer_contents
            !! Returns the contents of the solution buffer.
        procedure, public :: get_previous_error_norm => vsi_get_prev_err_norm
            !! Gets the norm of the previous step's error estimate.
        procedure, public :: set_previous_error_norm => vsi_set_prev_err_norm
            !! Sets the norm of the previous step's error estimate.
    end type

    interface
        subroutine variable_step_attempt(this, sys, h, x, y, yn, en, xprev, &
            yprev, fprev, err)
            !! Defines a routine meant to attempt a single integration step.
            use iso_fortran_env
            use ferror
            import variable_step_integrator
            import ode_container
            class(variable_step_integrator), intent(inout) :: this
                !! The variable_step_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the N
                !! dependent variables.
            real(real64), intent(out), dimension(:) :: yn
                !! An N-element array where the values of the dependent 
                !! variables at x + h will be written.
            real(real64), intent(out), dimension(:) :: en
                !! An N-element array where the error estimates for each 
                !! equation will be written.
            real(real64), intent(in), optional, dimension(:) :: xprev
                !! An optional M-element array containing the previous M values
                !! of the independent variable where M is the order of the 
                !! method.  This is typically only used for multi-step methods.
                !! In single-step methods, this parameter is typically not
                !! needed.
            real(real64), intent(in), optional, dimension(:,:) :: yprev
                !! An optional M-by-N matrix containing the previous M arrays of
                !! dependent variables, where M is the order of the method.  As
                !! with xprev, this parameter is typically used for multi-step
                !! methods.  In single-step methods, this parameter is typically
                !! not needed.
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
                !! An optional M-by-N matrix containing the previous M arrays of
                !! ODE (function) values.  As with xprev and yprev, M is the
                !! order of the method, and this parameter is typically used for
                !! multi-step methods.  In single-step methods, this parameter 
                !! is typically not needed.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to 
                !! provide error handling.
        end subroutine

        subroutine variable_step_action(this, x, xn, y, yn)
            !! Defines a routine for performing any actions upon completion of
            !! a successful step.
            use iso_fortran_env
            import variable_step_integrator
            class(variable_step_integrator), intent(inout) :: this
                !! The variable_step_integrator object.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in) :: xn
                !! The new value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the N
                !! dependent variables.
            real(real64), intent(in), dimension(:) :: yn
                !! An N-element array containing the new values of the N
                !! dependent variables.
        end subroutine

        subroutine variable_step_interpolation(this, xprev, yprev, xnew, x, y, &
            err)
            !! Defines a routine for providing interpolation services between
            !! integration points.
            use iso_fortran_env
            use ferror
            import variable_step_integrator
            class(variable_step_integrator), intent(in) :: this
                !! The variable_step_integrator object.
            real(real64), intent(in) :: xprev
                !! The previoius value of the independent variable.
            real(real64), intent(in), dimension(:) :: yprev
                !! An N-element array containing the previous values of the
                !! N dependent variables.
            real(real64), intent(in) :: xnew
                !! The new value of the independent variable.
            real(real64), intent(in) :: x
                !! The value at which to perform the interpolation.
            real(real64), intent(out), dimension(:) :: y
                !! An N-element array where the interpolated values will be
                !! written.
            class(errors), intent(inout), optional, target :: err
                !! An optional errors-based object that if provided 
                !! can be used to retrieve information relating to any errors 
                !! encountered during execution. If not provided, a default 
                !! implementation of the errors class is used internally to 
                !! provide error handling.
        end subroutine

        function next_step_size_calculator(this, hn, en, enm1) result(rst)
            !! Defines a routine for computing the next step size to attempt.
            use iso_fortran_env
            import variable_step_integrator
            class(variable_step_integrator), intent(inout) :: this
                !! The variable_step_integrator object.
            real(real64), intent(in) :: hn
                !! The current step size.
            real(real64), intent(in) :: en
                !! The norm of the error for the current step.
            real(real64), intent(in) :: enm1
                !! The norm of the error from the previous step.
            real(real64) :: rst
                !! The next step size to try.
        end function
    end interface

contains
! ------------------------------------------------------------------------------
pure function vsi_get_safety_factor(this) result(rst)
    !! Gets a safety factor used to limit the predicted step size.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The safety factor value.
    rst = this%m_safetyfactor
end function

! --------------------
subroutine vsi_set_safety_factor(this, x)
    !! Sets a safety factor used to limit the predicted step size.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The safety factor value.
    this%m_safetyfactor = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_max_step(this) result(rst)
    !! Gets the maximum allowed step size.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The maximum step size.
    rst = this%m_maxstep
end function

! --------------------
subroutine vsi_set_max_step(this, x)
    !! Sets the maximum step size.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The maximum step size.
    this%m_maxstep = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_min_step(this) result(rst)
    !! Gets the minimum allowed step size.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The minimum step size.
    rst = this%m_minstep
end function

! --------------------
subroutine vsi_set_min_step(this, x)
    !! Sets the minimum allowed step size.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The minimum step size.
    this%m_minstep = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_max_iter_count(this) result(rst)
    !! Gets the maximum number of iterations per step allowed.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    integer(int32) :: rst
        !! The maximum iteration count.
    rst = this%m_maxitercount
end function

! --------------------
subroutine vsi_set_max_iter_count(this, x)
    !! Sets the maximum number of iterations per step allowed.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    integer(int32), intent(in) :: x
        !! The maximum iteration count.
    this%m_maxitercount = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_max_step_count(this) result(rst)
    !! Gets the maximum number of integration steps allowed.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    integer(int32) :: rst
        !! The maximum number of integration steps.
    rst = this%m_maxstepcount
end function

! --------------------
subroutine vsi_set_max_step_count(this, x)
    !! Sets the maximum number of integration steps allowed.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    integer(int32), intent(in) :: x
        !! The maximum number of integration steps.
    this%m_maxstepcount = x
end subroutine

! ------------------------------------------------------------------------------
subroutine vsi_append_to_buffer(this, x, y, err)
    !! Buffers a results set.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The independent variable value.
    real(real64), intent(in) :: y(:)
        !! An N-element array containing the solution values.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.  Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if @p y is not compatible with
        !!      the buffer size.

    ! Parameters
    integer(int32), parameter :: buffer = 1000

    ! Local Variables
    integer(int32) :: m, n, neqn, flag
    real(real64), allocatable, dimension(:,:) :: copy
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    n = neqn + 1

    ! Ensure the buffer is allocated
    if (.not.allocated(this%m_buffer)) then
        allocate(this%m_buffer(buffer, n), stat = flag)
        if (flag /= 0) then
            call report_memory_error(errmgr, "vsi_append_to_buffer", flag)
            return
        end if
        this%m_bufferCount = 0
    end if

    ! Push a value onto the end of the buffer
    if (size(this%m_buffer, 2) /= n) then
        call report_array_size_error(errmgr, "vsi_append_to_buffer", "buffer", &
            n, size(this%m_buffer, 2))
        return
    end if

    this%m_bufferCount = this%m_bufferCount + 1
    m = size(this%m_buffer, 1)
    if (this%m_bufferCount > m) then
        ! Reallocate the buffer
        allocate(copy(m, n), stat = flag, source = this%m_buffer)
        if (flag /= 0) then
            call report_memory_error(errmgr, "vsi_append_to_buffer", flag)
            return
        end if
        deallocate(this%m_buffer)
        allocate(this%m_buffer(m + buffer, n), stat = flag)
        this%m_buffer(1:m,:) = copy
    end if

    this%m_buffer(this%m_bufferCount, 1) = x
    this%m_buffer(this%m_bufferCount, 2:) = y

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_buffer_count(this) result(rst)
    !! Gets the number of entries into the solution buffer.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    integer(int32) :: rst
        !! The number of buffer entries.
    rst = this%m_bufferCount
end function

! ------------------------------------------------------------------------------
subroutine vsi_clear_buffer(this)
    !! Clears the results buffer.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    if (allocated(this%m_buffer)) deallocate(this%m_buffer)
    this%m_bufferCount = 0
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_step_size(this) result(rst)
    !! Gets the current step size.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The step size.
    rst = this%m_stepSize
end function

! --------------------
subroutine vsi_set_step_size(this, x)
    !! Sets the current step size.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The step size.
    this%m_stepSize = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_next_step_size(this) result(rst)
    !! Gets the next step size.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The step size.
    rst = this%m_nextStep
end function

! --------------------
subroutine vsi_set_next_step_size(this, x)
    !! Sets the next step size.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The step size.
    this%m_nextStep = x
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_get_respect_xmax(this) result(rst)
    !! Gets a value determining if the integrator should respect a
    !! hard limit in the independent variable range.  If false, the 
    !! integrator may step pass the limit.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    logical :: rst
        !! True if the integrator should respect the limiting value of
        !!  the independent variable; else, false.
    rst = this%m_respectXMax
end function

! --------------------
subroutine vsi_set_respect_xmax(this, x)
    !! Sets a value determining if the integrator should respect a
    !! hard limit in the independent variable range.  If false, the 
    !! integrator may step pass the limit.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    logical, intent(in) :: x
        !! True if the integrator should respect the limiting 
        !!  value of the independent variable; else, false.
    this%m_respectXMax = x
end subroutine

! ------------------------------------------------------------------------------
subroutine vsi_alloc_workspace(this, neqn, err)
    !! Initializes the integrator.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.

    ! Local Variables
    integer(int32) :: flag
    real(real64) :: default_rtol, default_atol
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Process
    default_atol = this%get_default_absolute_tolerance()
    default_rtol = this%get_default_relative_tolerance()
    if (allocated(this%m_ework)) then
        if (size(this%m_ework) /= neqn) then
            deallocate(this%m_ework)
            allocate(this%m_ework(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_ework(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! Ensure tolerance arrays are initialized appropriately as well
    if (allocated(this%m_rtol)) then
        if (size(this%m_rtol) /= neqn) then
            deallocate(this%m_rtol)
            allocate(this%m_rtol(neqn), stat = flag, source = default_rtol)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_rtol(neqn), stat = flag, source = default_rtol)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_atol)) then
        if (size(this%m_atol) /= neqn) then
            deallocate(this%m_atol)
            allocate(this%m_atol(neqn), stat = flag, source = default_atol)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_atol(neqn), stat = flag, source = default_atol)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "vsi_alloc_workspace", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine vsi_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, err)
    !! Takes one integration step.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in) :: xmax
        !! The upper integration limit.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the dependent 
        !! variables.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where the values of the dependent variables at 
        !! x + h will be written.
    real(real64), intent(in), optional, dimension(:) :: xprev
        !! An optional M-element array containing the previous M values
        !! of the independent variable where M is the order of the 
        !! method.  This is typically only used for multi-step methods.
        !! In single-step methods, this parameter is typically not
        !! needed.
    real(real64), intent(in), optional, dimension(:,:) :: yprev
        !! An optional M-by-N matrix containing the previous M arrays of
        !! dependent variables, where M is the order of the method.  As
        !! with xprev, this parameter is typically used for multi-step
        !! methods.  In single-step methods, this parameter is typically
        !! not needed.
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
        !! An optional M-by-N matrix containing the previous M arrays of
        !! ODE (function) values.  As with xprev and yprev, M is the
        !! order of the method, and this parameter is typically used for
        !! multi-step methods.  In single-step methods, this parameter 
        !! is typically not needed.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h, enorm, et
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    h = this%get_next_step_size()

    ! Input Checking
    if (size(yn) /= neqn) then
        call report_array_size_error(errmgr, "vsi_step", "yn", neqn, size(yn))
        return
    end if

    ! Ensure the proper workspaces are allocated
    if (.not.allocated(this%m_ework)) then
        call this%initialize(neqn, errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! Process
    i = 0
    do
        ! Attempt a step
        call this%attempt_step(sys, h, x, y, yn, this%m_ework, xprev, &
            yprev, fprev, errmgr)
        if (errmgr%has_error_occurred()) return
        
        ! Compute the normalized error
        enorm = norm2( &
            this%m_ework / (neqn * (this%m_atol + &
                max(maxval(abs(y)), maxval(abs(yn))) * this%m_rtol)) &
        )
        if (enorm <= 1.0d0) call this%set_step_size(h)
            
        ! Compute a new step size
        h = this%compute_next_step_size(h, enorm, this%get_previous_error_norm())

        ! Check to see if the step size is too small
        if (abs(h) < abs(this%get_min_step_size())) then
            call report_step_size_too_small(errmgr, "vsi_step", x, h)
            return
        end if

        ! Do we need to limit the step size to not overstep a limiting x value?
        if (this%get_respect_x_max() .and. &
            abs(x + h) > abs(xmax)) &
        then
            h = xmax - x
        end if
        
        ! Is the step successful (enorm is normalized to the error tolerances
        ! such that a value less or equal to 1 is successful; else, keep going)
        if (enorm <= 1.0d0) exit

        ! Update the iteration counter
        i = i + 1
        if (i > this%get_max_per_step_iteration_count()) then
            call report_excessive_iterations(errmgr, "vsi_step", i, x)
            return
        end if
    end do

    ! Store error values from this step
    call this%set_previous_error_norm(enorm)

    ! Store the updated step size
    call this%set_next_step_size(h)

    ! Perform any actions needed on a successful step
    call this%on_successful_step(x, x + this%get_step_size(), y, yn)

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
pure function vsi_estimate_first_step(this, xo, xf, yo, fo) &
    result(rst)
    !! Computes an estimate to the first step size based upon the initial 
    !! function values.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: xo
        !! The initial value of the independent variable.
    real(real64), intent(in) :: xf
        !! The final value of the independent variable.
    real(real64), intent(in), dimension(:) :: yo
        !! An N-element array containing the initial values.
    real(real64), intent(in), dimension(:) :: fo
        !! An N-element array containing the initial function values.
    real(real64) :: rst
        !! An estimate on the initial step size.

    ! Local Variables
    real(real64) :: h1, h2
    
    ! Process
    h1 = 0.5d0 * (xf - xo)
    h2 = this%get_max_step_size()
    rst = sign(min(abs(h1), abs(h2)), h1)
end function

! ------------------------------------------------------------------------------
pure function vsi_get_default_rel_tol(this) result(rst)
    !! Gets the default relative error tolerance.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The tolerance value.
    rst = 1.0d-6
end function

! ------------------------------------------------------------------------------
pure function vsi_get_default_abs_tol(this) result(rst)
    !! Gets the default absolute error tolerance.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The tolerance value.
    rst = 1.0d-6
end function

! ------------------------------------------------------------------------------
function vsi_get_buffer_contents(this) result(rst)
    !! Returns the contents of the solution buffer.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The buffer contents.

    ! Local Variables
    integer(int32) :: m, n
    
    ! Initialization
    m = this%get_buffer_size()
    n = size(this%m_buffer, 2)

    ! Process
    if (allocated(this%m_buffer)) then
        allocate(rst(m, n), source = this%m_buffer(:m,:))
    else
        allocate(rst(0, 0))
    end if
end function

! ------------------------------------------------------------------------------
pure function vsi_get_prev_err_norm(this) result(rst)
    !! Gets the norm of the previous step's error estimate.
    class(variable_step_integrator), intent(in) :: this
        !! The variable_step_integrator object.
    real(real64) :: rst
        !! The error norm.
    rst = this%m_enormPrev
end function

! --------------------
subroutine vsi_set_prev_err_norm(this, x)
    !! Sets the norm of the previous step's error estimate.
    class(variable_step_integrator), intent(inout) :: this
        !! The variable_step_integrator object.
    real(real64), intent(in) :: x
        !! The error norm.
    this%m_enormPrev = x
end subroutine

! ------------------------------------------------------------------------------
end module