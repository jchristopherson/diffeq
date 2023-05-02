submodule (diffeq) diffeq_multistep_fixed
contains
! ------------------------------------------------------------------------------
module subroutine fms_alloc_workspace(this, neqn, err)
    ! Arguments
    class(fixed_multistep_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: n, flag
    character(len = :), allocatable :: errmsg

    ! Process
    n = this%get_order()
    if (allocated(this%m_buffer)) then
        if (size(this%m_buffer, 1) /= neqn .or. &
            size(this%m_buffer, 2) /= n) &
        then
            deallocate(this%m_buffer)
            allocate(this%m_buffer(neqn, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_buffer(neqn, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("fms_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module function fms_solver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(fixed_multistep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) ::i, j, j1, j2, npts, neqn, order, mn, flag
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

    call this%allocate_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Use a 4th order Runge-Kutta integrator to step into the problem
    rst(1,1) = x(1)
    rst(1,2:) = iv
    call sys%fcn(x(1), iv, this%m_buffer(:,order))
    j2 = order - 1
    do i = 2, mn
        j = i - 1
        h = x(j) - x(i)
        call starter%step(sys, h, x(j), rst(j,2:), rst(i,2:))
        call sys%fcn(x(i), rst(i,2:), this%m_buffer(:,j2))
        j2 = j2 - 1
    end do

    ! Finish the problem with the multi-step method
    j1 = 1
    j2 = order
    do i = order + 1, npts
        h = x(i) - x(j2)
        rst(i,1) = x(i)
        call this%step(sys, h, rst(j2,1), rst(j2,2:), rst(i,2:), &
            rst(j1:j2,1), rst(j1:j2,2:), this%m_buffer, err = errmgr)
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