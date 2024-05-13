module diffeq_runge_kutta
    use iso_fortran_env
    use diffeq_errors
    use diffeq_base
    implicit none
    private
    public :: runge_kutta_45
    public :: runge_kutta_23

    type, extends(single_step_integrator) :: runge_kutta_45
        !! The Dormand-Prince, Runge-Kutta integrator (4th/5th order).
        real(real64), private, allocatable, dimension(:) :: k2
        real(real64), private, allocatable, dimension(:) :: k3
        real(real64), private, allocatable, dimension(:) :: k4
        real(real64), private, allocatable, dimension(:) :: k5
        real(real64), private, allocatable, dimension(:) :: k6
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
        real(real64), private, allocatable, dimension(:) :: rc4
        real(real64), private, allocatable, dimension(:) :: rc5
    contains
        procedure, public :: get_order => rk45_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk45_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, private :: initialize => rk45_init
            !! Initializes private storage arrays for the integrator.
        procedure, private :: initialize_interp => rk45_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: attempt_step => rk45_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rk45_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rk45_interp
            !! Performs the interpolation.
    end type

    type, extends(single_step_integrator) :: runge_kutta_23
        !! The Bogacki-Shampine integrator (3rd/2nd order).
        real(real64), private, allocatable, dimension(:) :: k2
        real(real64), private, allocatable, dimension(:) :: k3
        real(real64), private, allocatable, dimension(:) :: k4
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
    contains
        procedure, private :: initialize => rk32_init
            !! Initializes private storage arrays for the integrator.
        procedure, private :: initialize_interp => rk32_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: get_order => rk32_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk32_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: attempt_step => rk32_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rk32_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rk32_interp
            !! Performs the interpolation.
    end type

contains
! ******************************************************************************
! RUNGE_KUTTA_45
! ------------------------------------------------------------------------------
pure function rk45_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    integer(int32) :: rst
        !! The order.
    rst = 5
end function

! ------------------------------------------------------------------------------
pure function rk45_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine rk45_init(this, neqn)
    !! Initializes private storage arrays for the integrator.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%k2)) then
        if (size(this%k2) == neqn) then
            ! All is good - go ahead and return
            return
        else
            ! The arrays are not sized properly
            deallocate(this%k2)
            deallocate(this%k3)
            deallocate(this%k4)
            deallocate(this%k5)
            deallocate(this%k6)
        end if
    end if

    ! If we're here, we need to allocate the storage arrays
    allocate( &
        this%k2(neqn), &
        this%k3(neqn), &
        this%k4(neqn), &
        this%k5(neqn), &
        this%k6(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_init_interp(this, neqn)
    !! Allocates storage for the interpolation process.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%rc1)) then
        if (size(this%rc1) == neqn) then
            ! All is good
            return
        else
            deallocate(this%rc1)
            deallocate(this%rc2)
            deallocate(this%rc3)
            deallocate(this%rc4)
            deallocate(this%rc5)
        end if
    end if

    allocate( &
        this%rc1(neqn), &
        this%rc2(neqn), &
        this%rc3(neqn), &
        this%rc4(neqn), &
        this%rc5(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_attempt_step(this, sys, h, x, y, f, yn, fn, yerr)
    use diffeq_dprk45_constants
    !! Attempts an integration step for this integrator.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives
        !! at x.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where this routine will write the next
        !! solution estimate at x + h.
    real(real64), intent(out), dimension(:) :: fn
        !! An N-element array where this routine will write the next
        !! derivative estimate at x + h.
    real(real64), intent(out), dimension(:) :: yerr
        !! An N-element array where this routine will write an estimate
        !! of the error in each equation.

    ! Local Variables
    integer(int32) :: n

    ! Initialization
    n = size(y)
    call this%initialize(n)

    ! Process
    ! k1 = f as the derivatives were computed from the previous step

    yn = y + h * a21
    call sys%fcn(x + h * c2, yn, this%k2)

    yn = y + h * (a31 * f + a32 * this%k2)
    call sys%fcn(x + h * c3, yn, this%k3)

    yn = y + h * (a41 * f + a42 * this%k2 + a43 * this%k3)
    call sys%fcn(x + h * c4, yn, this%k4)

    yn = y + h * (a51 * f + a52 * this%k2 + a53 * this%k3 + a54 * this%k4)
    call sys%fcn(x + h * c5, yn, this%k5)

    yn = y + h * (a61 * f + a62 * this%k2 + a63 * this%k3 + a64 * this%k4 + &
        a65 * this%k5)
    call sys%fcn(x + h * c6, yn, this%k6)

    yn = y + h * (a71 * f + a73 * this%k3 + a74 * this%k4 + a75 * this%k5 + a76 * this%k6)
    call sys%fcn(x + h * c7, yn, fn)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * this%k2 + e3 * this%k3 + e4 * this%k4 + &
        e5 * this%k5 + e6 * this%k6 + e7 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_set_up_interp(this, dense, x, xn, y, yn, f, fn)
    use diffeq_dprk45_constants
    !! Sets up the interpolation process.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    logical, intent(in) :: dense
        !! Determines if dense output is requested (true); else, false.
    real(real64), intent(in) :: x
        !! The previous value of the independent variable.
    real(real64), intent(in) :: xn
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution at x.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the derivatives at x.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: h, ydiff, bspl

    ! Initialization
    n = size(y)
    call this%initialize_interp(n)
    h = xn - x

    ! Process
    do i = 1, n
        this%rc1(i) = y(i)
        ydiff = yn(i) - y(i)
        this%rc2(i) = ydiff
        bspl = h * f(i) - ydiff
        this%rc3(i) = bspl
        this%rc4(i) = ydiff - h * fn(i) - bspl
        this%rc5 = h * (d1 * f(i) + d3 * this%k3(i) + d4 * this%k4(i) + &
            d5 * this%k5(i) + d6 * this%k6(i) + d7 * fn(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which to compute
        !! the interpolation.
    real(real64), intent(in) :: xn
        !! The previous value of the independent variable at which the
        !! solution is computed.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.
    real(real64), intent(in) :: xn1
        !! The value of the independent variable at xn + h.
    real(real64), intent(in), dimension(:) :: yn1
        !! An N-element array containing the solution at xn + h.
    real(real64), intent(in), dimension(:) :: fn1
        !! An N-element array containing the derivatives at xn + h.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where this routine will write the 
        !! solution values interpolated at x.

    ! Local Variables
    real(real64) :: h, s, s1

    ! Process
    h = xn1 - xn
    s = (x - xn) / h
    s1 = 1.0d0 - s
    y = this%rc1 + s * (this%rc2 + s1 * (this%rc3 + &
        s * (this%rc4 + s1 * this%rc5)))
end subroutine

! ******************************************************************************
! RUNGE_KUTTA_23
! ------------------------------------------------------------------------------
pure function rk32_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    integer(int32) :: rst
        !! The order.
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function rk32_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine rk32_init(this, neqn)
    !! Initializes private storage arrays for the integrator.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%k2)) then
        if (size(this%k2) == neqn) then
            ! All is good - go ahead and return
            return
        else
            ! The arrays are not sized properly
            deallocate(this%k2)
            deallocate(this%k3)
            deallocate(this%k4)
        end if
    end if

    ! If we're here, we need to allocate the storage arrays
    allocate( &
        this%k2(neqn), &
        this%k3(neqn), &
        this%k4(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk32_init_interp(this, neqn)
    !! Allocates storage for the interpolation process.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%rc1)) then
        if (size(this%rc1) == neqn) then
            ! All is good
            return
        else
            deallocate(this%rc1)
            deallocate(this%rc2)
            deallocate(this%rc3)
        end if
    end if

    allocate( &
        this%rc1(neqn), &
        this%rc2(neqn), &
        this%rc3(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk32_attempt_step(this, sys, h, x, y, f, yn, fn, yerr)
    use diffeq_bsrk32_constants
    !! Attempts an integration step for this integrator.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives
        !! at x.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where this routine will write the next
        !! solution estimate at x + h.
    real(real64), intent(out), dimension(:) :: fn
        !! An N-element array where this routine will write the next
        !! derivative estimate at x + h.
    real(real64), intent(out), dimension(:) :: yerr
        !! An N-element array where this routine will write an estimate
        !! of the error in each equation.

    ! Local Variables
    integer(int32) :: n

    ! Initialization
    n = size(y)
    call this%initialize(n)

    ! Process
    ! k1 = f as the derivatives were computed from the previous step

    yn = y + h * a21
    call sys%fcn(x + h * c2, yn, this%k2)

    yn = y + h * (a31 * f + a32 * this%k2)
    call sys%fcn(x + h * c3, yn, this%k3)

    yn = y + h * (a41 * f + a42 * this%k2 + a43 * this%k3)
    call sys%fcn(x + h * c4, yn, fn)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * this%k2 + e3 * this%k3 + e4 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk32_set_up_interp(this, dense, x, xn, y, yn, f, fn)
    use diffeq_bsrk32_constants
    !! Sets up the interpolation process.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    logical, intent(in) :: dense
        !! Determines if dense output is requested (true); else, false.
    real(real64), intent(in) :: x
        !! The previous value of the independent variable.
    real(real64), intent(in) :: xn
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution at x.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the derivatives at x.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.

    ! Local Variables
    integer(int32) :: n
    real(real64) :: h

    ! Initialization
    h = xn - x
    n = size(y)
    call this%initialize_interp(n)

    ! Process
    this%rc1 = -(y - yn + f * h) / h**2
    this%rc2 = f - 2.0d0 * x * this%rc1
    this%rc3 = y - this%rc1 * x**2 - this%rc2 * x
end subroutine
! ------------------------------------------------------------------------------
subroutine rk32_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which to compute
        !! the interpolation.
    real(real64), intent(in) :: xn
        !! The previous value of the independent variable at which the
        !! solution is computed.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.
    real(real64), intent(in) :: xn1
        !! The value of the independent variable at xn + h.
    real(real64), intent(in), dimension(:) :: yn1
        !! An N-element array containing the solution at xn + h.
    real(real64), intent(in), dimension(:) :: fn1
        !! An N-element array containing the derivatives at xn + h.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where this routine will write the 
        !! solution values interpolated at x.

    ! Process
    y = this%rc1 * x**2 + this%rc2 * x + this%rc3
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module