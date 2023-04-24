submodule (diffeq) diffeq_fs_integrator
contains
! ------------------------------------------------------------------------------
module function fsi_solver(this, x, iv, err) result(rst)
    ! Arguments
    class(fixed_step_integrator), intent(in) :: this
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: npts, neqn, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    neqn = size(iv)

    ! Input Checking
    if (.not.associated(this%fcn)) go to 20
    if (npts < 2) go to 30


    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    rst(1,1) = x(1)
    rst(1,2:) = iv
    do i = 2, npts
        ! Take the step
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    return

    ! No Function Defined Error
20  continue
    return

    ! Independent Variable Array Size Error
30  continue
    return
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule