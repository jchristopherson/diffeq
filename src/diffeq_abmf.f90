submodule (diffeq) diffeq_abmf
contains
! ------------------------------------------------------------------------------
pure module function afi_get_order(this) result(rst)
    class(adams_fixed_integerator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
module subroutine afi_alloc_workspace(this, neqn, err)
    ! Arguments
    class(adams_fixed_integerator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: n, flag
    character(len = :), allocatable :: errmsg

    ! Process
    n = this%get_order()
    this%m_first = .true.
    if (allocated(this%m_work)) then
        if (size(this%m_work, 1) /= neqn .or. size(this%m_work, 2) /= n) then
            allocate(this%m_work(neqn, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_work(neqn, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("afi_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine afi_shift_workspace(this)
    ! Arguments
    class(adams_fixed_integerator), intent(inout) :: this

    ! Local Variables
    integer(int32) :: j, n

    ! Shift each column over by 1 allowing the last column to fall off
    n = size(this%m_work, 2)
    do j = n, 2, -1
        this%m_work(:,j) = this%m_work(:,j-1)
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine afi_step(this, sys, h, x, y, yn, xprev, yprev, err)
    ! Arguments
    class(adams_fixed_integerator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    class(errors), intent(inout), optional, target :: err

    ! Model Constants
    real(real64), parameter :: a1 = 55.0d0 / 24.0d0
    real(real64), parameter :: a2 = -59.0d0 / 24.0d0
    real(real64), parameter :: a3 = 37.0d0 / 24.0d0
    real(real64), parameter :: a4 = -9.0d0 / 24.0d0
    real(real64), parameter :: b1 = 9.0d0 / 24.0d0
    real(real64), parameter :: b2 = 19.0d0 / 24.0d0
    real(real64), parameter :: b3 = -5.0d0 / 24.0d0
    real(real64), parameter :: b4 = 1.0d0 / 24.0d0

    ! Local Variables
    integer(int32) :: i, j, n
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
    n = this%get_order()

    ! Input Checking
    if (size(yn) /= neqn) go to 10
    if (.not.present(xprev) .or. .not.present(yprev)) go to 20
    if (size(xprev) < n) go to 30
    if (size(yprev, 1) < n .or. size(yprev, 2) /= neqn) go to 40

    ! Allocate the workspace
    call this%allocate_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Process
    if (this%m_first) then
        ! First time through - evaluate functions and store in the workspace
        j = n
        do i = 2, n
            call sys%fcn(xprev(j), yprev(j,:), this%m_work(:,i))
            j = j - 1
        end do
        this%m_first = .false.
    end if

    ! Compute the Adams-Bashforth predictor
    yn = y + h * ( &
        a1 * this%m_work(:,1) + &
        a2 * this%m_work(:,2) + &
        a3 * this%m_work(:,3) + &
        a4 * this%m_work(:,4) &
    )
    call this%shift()
    call sys%fcn(x + h, yn, this%m_work(:,1))

    ! Compute the Adams-Moulton corrector
    yn = y + h * ( &
        b1 * this%m_work(:,1) + &
        b2 * this%m_work(:,2) + &
        b3 * this%m_work(:,3) + &
        b4 * this%m_work(:,4) &
    )
    call sys%fcn(x + h, yn, this%m_work(:,1))

    ! End
    return

    ! YN array size error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to have ", neqn, &
        " elements, but was found to have ", size(yn), " elements."
    call errmgr%report_error("afi_step", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! XPREV or YPREV not input
20  continue
    return

    ! XPREV wrong size
30  continue
    return

    ! YPREV wrong size
40  continue
    return

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule