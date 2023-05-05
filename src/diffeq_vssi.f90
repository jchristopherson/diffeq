submodule (diffeq) diffeq_vssi
contains
! ------------------------------------------------------------------------------
module function vssi_solve(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(variable_singlestep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: nx
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nx = size(x)

    ! Input Checking
    if (nx < 2) go to 10
    if (.not.associated(sys%fcn)) go to 20
    if (abs(maxval(x) - minval(x)) < epsilon(1.0d0)) go to 30

    ! Process
    if (nx == 2) then
        rst = this%solve_driver(sys, x, iv, errmgr)
    else
        rst = this%dense_solve_driver(sys, x, iv, errmgr)
    end if
    if (errmgr%has_error_occurred()) return

    ! End
    return

    ! X isn't sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The independent variable array must have at " // &
        "least 2 elements, but was found to have ", nx, "."
    call errmgr%report_error("vssi_solve", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! No ODE is defined
20  continue
    call errmgr%report_error("vssi_solve", "No ODE routine defined.", &
        DIFFEQ_NULL_POINTER_ERROR)
    return

    ! XMAX - XMIN == 0 error
30  continue
    call errmgr%report_error("vssi_solve", "MAX(X) - MIN(X) is zero.", &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
module function vssi_solve_driver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(variable_singlestep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout) :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Parameters
    real(real64), parameter :: sml = 1.0d2 * epsilon(1.0d0)

    ! Local Variables
    integer(int32) :: i, neqn, flag, order, ns
    real(real64) :: xn, xmax
    real(real64), allocatable, dimension(:) :: yn, yn1
    character(len = :), allocatable :: errmsg

    ! Initialization
    xmax = x(2)
    xn = x(1)
    neqn = size(iv)
    order = this%get_order()

    ! Memory Allocation
    allocate(yn(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(yn1(neqn), stat = flag)
    if (flag /= 0) go to 10

    ! Store the initial conditions
    call this%buffer_results(x(1), iv, err)
    if (err%has_error_occurred()) return

    ! Provide an initial step size estimate
    call this%set_next_step_size(min(0.5d0 * (x(2) - x(1)), &
        this%get_max_step_size()))

    ! Cycle until complete
    i = 0
    do
        ! Take the step
        call this%step(sys, xn, xmax, yn, yn1, err = err)
        if (err%has_error_occurred()) return

        ! Update xn and yn
        xn = xn + this%get_step_size()
        yn = yn1

        ! Store the solution
        call this%buffer_results(xn, yn, err)
        if (err%has_error_occurred()) return

        ! Break if |XN| > |XMAX|
        if (abs(xn) >= abs(xmax) .or. abs(this%get_step_size()) < sml) exit

        ! Iteration Counter
        i = i + 1
        if (i > this%get_max_integration_step_count()) go to 20
    end do

    ! Output the results
    ns = this%get_buffer_size()
    rst = this%m_buffer(1:ns,:)

    ! End
    return

    ! Memory Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("vssi_solve_driver", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Too many steps
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The allowable number of integration steps " // &
        "was exceeded at x = ", xn, "."
    call err%report_error("vssi_solve_driver", trim(errmsg), &
        DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, E10.3, A)
end function

! ------------------------------------------------------------------------------
module function vssi_dense_solve_driver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(variable_singlestep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout) :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Parameters
    real(real64), parameter :: sml = 1.0d2 * epsilon(1.0d0)

    ! Local Variables
    integer(int32) :: i, j, npts, neqn, flag
    real(real64) :: xn, xmax
    real(real64), allocatable, dimension(:) :: yn, yn1
    character(len = :), allocatable :: errmsg
    
    ! Initialization

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag == 0) allocate(yn(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(yn1(neqn), stat = flag)
    if (flag /= 0) go to 10

    ! Store the initial conditions
    rst(1,1) = x(1)
    rst(1,2:) = iv

    ! Cycle until complete
    i = 0
    do
        ! Take a step
        call this%step(sys, xn, xmax, yn, yn1, err = err)
        if (err%has_error_occurred()) return

        ! Interpolate as needed to achieve any intermediary solution points

        ! Update xn and yn

        ! Break if |XN| > |XMAX|
        if (abs(xn) >= abs(xmax) .or. abs(this%get_step_size()) < sml) exit

        ! Iteration Counter
        i = i + 1
        if (i > max(this%get_max_integration_step_count(), npts)) go to 20
    end do

    ! End
    return

    ! Memory Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("vssi_solve_driver", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Too many steps
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The allowable number of integration steps " // &
        "was exceeded at x = ", xn, "."
    call err%report_error("vssi_solve_driver", trim(errmsg), &
        DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, E10.3, A)
end function

! ------------------------------------------------------------------------------
end submodule