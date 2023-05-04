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

    ! Process
    if (nx == 2) then
        rst = this%solve_driver(sys, x, iv, errmgr)
    else
    end if

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
    return

    ! XMAX - XMIN == 0 error
30  continue
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
    integer(int32) :: neqn, flag, order, ns
    real(real64) :: xn, xmax
    real(real64), allocatable, dimension(:) :: yn, yn1

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
    call this%set_next_step_size(0.5d0 * (x(2) - x(1)))

    ! Cycle until complete
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
    end do

    ! Output the results
    ns = this%get_buffer_size()
    rst = this%m_buffer(1:ns,:)

    ! End
    return

    ! Memory Error
10  continue
    return
end function

! ------------------------------------------------------------------------------
end submodule