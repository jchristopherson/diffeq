module diffeq_abmf
    !! This module provides a fixed-step Adams-type integrator.
    use iso_fortran_env
    use diffeq_multistep_fixed
    use diffeq_base
    use diffeq_errors
    use ferror
    implicit none
    private
    public :: adams_fixed_integerator

    type, extends(fixed_multistep_integrator) :: adams_fixed_integerator
        !! Defines a fixed-step, 4th-order, Adams-Bashforth-Moulton PECE 
        !! integrator.
    contains
        procedure, public :: get_order => afi_get_order
            !! Returns the order of the integrator.
        procedure, public :: step => afi_step
            !! Takes a single fixed-size integration step.
    end type

contains
! ------------------------------------------------------------------------------
pure module function afi_get_order(this) result(rst)
    !! Returns the order of the integrator.
    class(adams_fixed_integerator), intent(in) :: this
        !! The adams_fixed_integerator object.
    integer(int32) :: rst
        !! The order of the integrator.
    rst = 4
end function

! ------------------------------------------------------------------------------
subroutine shift(x)
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
    !! Takes a single fixed-size integration step.
    class(adams_fixed_integerator), intent(inout) :: this
        !! The adams_fixed_integerator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The size of the step to take.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the
        !! dependent variables.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where the values of the dependent 
        !! variables at x + h will be written.
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
        !!
        !! - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if any of the input arrays are not
        !!      sized appropriately.
        !!
        !! - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if any of the input matrices are
        !!      not sized appropriately.
        !!
        !! - DIFFEQ_MISSING_ARGUMENT_ERROR: Occurs if xprev, yprev, and/or
        !!      fprev are not provided.

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
    integer(int32) :: i, j, n, neqn
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
    if (size(yn) /= neqn) then
        call report_array_size_error(errmgr, "afi_step", "yn", neqn, size(yn))
        return
    end if
    if (.not.present(yprev)) then
        call report_missing_argument(errmgr, "afi_step", "xprev")
        return
    end if
    if (.not.present(yprev)) then
        call report_missing_argument(errmgr, "afi_step", "yprev")
        return
    end if
    if (.not.present(fprev)) then
        call report_missing_argument(errmgr, "afi_step", "fprev")
        return
    end if
    if (size(xprev) < n) then
        call report_array_size_error(errmgr, "afi_step", "xprev", n, &
            size(xprev))
        return
    end if
    if (size(yprev, 1) < n .or. size(yprev, 2) /= neqn) then
        call report_matrix_size_error(errmgr, "afi_step", "yprev", n, neqn, &
            size(yprev, 1), size(yprev, 2))
        return
    end if
    if (size(fprev, 1) < n .or. size(fprev, 2) /= neqn) then
        call report_matrix_size_error(errmgr, "afi_step", "fprev", n, neqn, &
            size(fprev, 1), size(fprev, 2))
        return
    end if

    ! Compute the Adams-Bashforth predictor
    yn = y + h * ( &
        a1 * fprev(:,1) + &
        a2 * fprev(:,2) + &
        a3 * fprev(:,3) + &
        a4 * fprev(:,4) &
    )
    call shift(fprev)
    call sys%ode(x + h, yn, fprev(:,1))

    ! Compute the Adams-Moulton corrector
    yn = y + h * ( &
        b1 * fprev(:,1) + &
        b2 * fprev(:,2) + &
        b3 * fprev(:,3) + &
        b4 * fprev(:,4) &
    )
    call sys%ode(x + h, yn, fprev(:,1))

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
end module