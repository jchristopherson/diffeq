submodule (diffeq) diffeq_multistep_fixed
contains
! ------------------------------------------------------------------------------
module function fms_solver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(fixed_multistep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) ::i, j1, j2, npts, neqn, order, mn, flag
    real(real64) :: h
    type(rk4_fixed_integrator) :: starter
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    neqn = size(iv)
    order = this%get_order()
    mn = min(order, npts)

    ! Input Checking
    if (.not.associated(sys%fcn)) go to 20
    if (npts < 2) go to 30

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag /= 0) go to 10

    ! Use a 4th order Runge-Kutta integrator to step into the problem
    rst(1:mn,:) = starter%solve(sys, x(1:mn), iv, errmgr)

    ! Finish the problem with the multi-step method
    j1 = 1
    j2 = order
    do i = order + 1, npts
        h = x(i) - x(j2)
        rst(i,1) = x(i)
        call this%step(sys, h, rst(j2,1), rst(j2,2:), rst(i,2:), &
            rst(j1:j2,1), rst(j1:j2,2:), errmgr)
        if (errmgr%has_error_occurred()) return
        j1 = j1 + 1
        j2 = j2 + 1
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("fms_solver", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! No Function Defined Error
20  continue
    call errmgr%report_error("fms_solver", "The ODE routine is not defined.", &
        DIFFEQ_NULL_POINTER_ERROR)
    return

    ! Independent Variable Array Size Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) &
        "There must be at least 2 solution points defined, but ", npts, &
        " were found."
    call errmgr%report_error("fms_solver", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
end submodule