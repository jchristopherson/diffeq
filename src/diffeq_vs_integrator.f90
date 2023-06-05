submodule (diffeq) diffeq_vs_integrator
contains
! ------------------------------------------------------------------------------
pure module function vsi_get_safety_factor(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_safetyfactor
end function

! --------------------
module subroutine vsi_set_safety_factor(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_safetyfactor = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_max_step(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_maxstep
end function

! --------------------
module subroutine vsi_set_max_step(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_maxstep = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_min_step(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_minstep
end function

! --------------------
module subroutine vsi_set_min_step(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_minstep = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_max_iter_count(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_maxitercount
end function

! --------------------
module subroutine vsi_set_max_iter_count(this, x)
    class(variable_step_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_maxitercount = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_max_step_count(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_maxstepcount
end function

! --------------------
module subroutine vsi_set_max_step_count(this, x)
    class(variable_step_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_maxstepcount = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine vsi_append_to_buffer(this, x, y, err)
    ! Arguments
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    real(real64), intent(in) :: y(:)
    class(errors), intent(inout), optional, target :: err

    ! Parameters
    integer(int32), parameter :: buffer = 1000

    ! Local Variables
    integer(int32) :: m, n, neqn, flag
    real(real64), allocatable, dimension(:,:) :: copy
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
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
        if (flag /= 0) go to 10
        this%m_bufferCount = 0
    end if

    ! Push a value onto the end of the buffer
    if (size(this%m_buffer, 2) /= n) go to 20

    this%m_bufferCount = this%m_bufferCount + 1
    m = size(this%m_buffer, 1)
    if (this%m_bufferCount > m) then
        ! Reallocate the buffer
        allocate(copy(m, n), stat = flag, source = this%m_buffer)
        if (flag /= 0) go to 10
        deallocate(this%m_buffer)
        allocate(this%m_buffer(m + buffer, n), stat = flag)
        this%m_buffer(1:m,:) = copy
    end if

    this%m_buffer(this%m_bufferCount, 1) = x
    this%m_buffer(this%m_bufferCount, 2:) = y

    ! End
    return

    ! Memory Error Issue
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("vsi_append_to_buffer", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Array Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The array size (", size(y), &
        ") is incompatible with the buffer (", size(this%m_buffer, 2) - 1, ")."
    call errmgr%report_error("vsi_append_to_buffer", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_buffer_count(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_bufferCount
end function

! ------------------------------------------------------------------------------
module subroutine vsi_clear_buffer(this)
    class(variable_step_integrator), intent(inout) :: this
    if (allocated(this%m_buffer)) deallocate(this%m_buffer)
    this%m_bufferCount = 0
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_step_size(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_stepSize
end function

! --------------------
module subroutine vsi_set_step_size(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_stepSize = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_next_step_size(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_nextStep
end function

! --------------------
module subroutine vsi_set_next_step_size(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_nextStep = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_respect_xmax(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    logical :: rst
    rst = this%m_respectXMax
end function

! --------------------
module subroutine vsi_set_respect_xmax(this, x)
    class(variable_step_integrator), intent(inout) :: this
    logical, intent(in) :: x
    this%m_respectXMax = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine vsi_alloc_workspace(this, neqn, err)
    ! Arguments
    class(variable_step_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: flag
    real(real64) :: default_rtol, default_atol
    character(len = :), allocatable :: errmsg
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
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("vsi_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine vsi_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(variable_step_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, xmax
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h, enorm, et
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    h = this%get_next_step_size()

    ! Input Checking
    if (size(yn) /= neqn) go to 10

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
        h = this%compute_next_step_size(h, enorm, this%m_enormPrev)

        ! Check to see if the step size is too small
        if (abs(h) < abs(this%get_min_step_size())) go to 20

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
        if (i > this%get_max_per_step_iteration_count()) go to 30
    end do

    ! Store error values from this step
    this%m_enormPrev = enorm

    ! Store the updated step size
    call this%set_next_step_size(h)

    ! Perform any actions needed on a successful step
    call this%on_successful_step(x, x + this%get_step_size(), y, yn)

    ! End
    return

    ! YN array size error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to have ", neqn, &
        " elements, but was found to have ", size(yn), " elements."
    call errmgr%report_error("vsi_step", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Step size too small error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "A step size of ", h, " was encountered at x = ", x, &
        ", and is too small to continue integration."
    call errmgr%report_error("vsi_step", trim(errmsg), &
        DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR)
    return

    ! Iteration count exceeded error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 102) "The allowable iteration count was exceeded at x = ", &
        x, "."
    call errmgr%report_error("vsi_step", trim(errmsg), &
        DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, E10.3, A, E10.3, A)
102 format(A, E10.3, A)
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_estimate_first_step(this, xo, xf, yo, fo) &
    result(rst)
    ! Arguments
    class(variable_step_integrator), intent(in) :: this
    real(real64), intent(in) :: xo, xf
    real(real64), intent(in), dimension(:) :: yo, fo
    real(real64) :: rst

    ! Local Variables
    real(real64) :: h1, h2
    
    ! Process
    h1 = 0.5d0 * (xf - xo)
    h2 = this%get_max_step_size()
    rst = min(h1, h2)
end function

! ------------------------------------------------------------------------------
pure module function vsi_get_default_rel_tol(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = 1.0d-6
end function

! ------------------------------------------------------------------------------
pure module function vsi_get_default_abs_tol(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = 1.0d-6
end function

! ------------------------------------------------------------------------------
end submodule