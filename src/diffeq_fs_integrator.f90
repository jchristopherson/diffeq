submodule (diffeq) diffeq_fs_integrator
contains
! ------------------------------------------------------------------------------
module function fsi_solver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(fixed_step_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: npts, neqn, flag
    real(real64) :: h
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

    ! Input Checking
    if (.not.associated(sys%fcn)) go to 20
    if (npts < 2) go to 30

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    rst(1,1) = x(1)
    rst(1,2:) = iv
    do i = 2, npts
        ! Compute the solution
        j = i - 1
        h = x(i) - x(j)
        call this%step(sys, h, rst(j,1), rst(j,2:), rst(i,2:))
        rst(i,1) = x(i)
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("fsi_solver", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! No Function Defined Error
20  continue
    call errmgr%report_error("fsi_solver", "The ODE routine is not defined.", &
        DIFFEQ_NULL_POINTER_ERROR)
    return

    ! Independent Variable Array Size Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) &
        "There must be at least 2 solution points defined, but ", npts, &
        " were found."
    call errmgr%report_error("fsi_solver", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule