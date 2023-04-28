submodule (diffeq) diffeq_abmf
contains
! ------------------------------------------------------------------------------
pure module function afi_get_order(this) result(rst)
    class(adams_fixed_integerator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
module subroutine shift(x)
    ! Arguments
    real(real64), intent(inout), dimension(:,:) :: x

    ! Local Variables
    integer(int32) :: j, n

    ! Shift each column over by 1 allowing the last column to fall off
    n = size(x, 2)
    do j = n, 2, -1
        x(:,j) = x(:,j-1)
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine afi_step(this, sys, h, x, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(adams_fixed_integerator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
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

    ! Compute the Adams-Bashforth predictor
    yn = y + h * ( &
        a1 * fprev(:,1) + &
        a2 * fprev(:,2) + &
        a3 * fprev(:,3) + &
        a4 * fprev(:,4) &
    )
    call shift(fprev)
    call sys%fcn(x + h, yn, fprev(:,1))

    ! Compute the Adams-Moulton corrector
    yn = y + h * ( &
        b1 * fprev(:,1) + &
        b2 * fprev(:,2) + &
        b3 * fprev(:,3) + &
        b4 * fprev(:,4) &
    )
    call sys%fcn(x + h, yn, fprev(:,1))

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

    ! XPREV, YPREV, or FPREV not input
20  continue
    call errmgr%report_error("afi_step", &
        "All optional array arguments must be supplied.", &
        DIFFEQ_MISSING_ARGUMENT_ERROR)
    return

    ! XPREV wrong size
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The independent variable history array was " // &
        "expected to have ", n, " elements, but was found to have ", &
        size(xprev), "."
    call errmgr%report_error("afi_step", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! YPREV wrong size
40  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The dependent variable history matrix was " // &
        "expected to be (", neqn, "-by-", n, "), but was found to be (", &
        size(yprev, 1), "-by-", size(yprev, 2), ")."
    call errmgr%report_error("afi_step", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! FPREV wrong size
50  continue
    write(errmsg, 101) "The ODE result history matrix was " // &
        "expected to be (", neqn, "-by-", n, "), but was found to be (", &
        size(yprev, 1), "-by-", size(yprev, 2), ")."
    call errmgr%report_error("afi_step", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
101 format(A, I0, A, I0, A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule