module diffeq_runge_kutta
    use iso_fortran_env
    use diffeq_errors
    use diffeq_base
    implicit none
    private
    public :: runge_kutta_45
    public :: runge_kutta_23
    public :: runge_kutta_853
    public :: diagonally_implicit_runge_kutta

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
        procedure, public :: pre_step_action => rk45_pre_step
            !! Performs any pre-step actions.
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
        procedure, public :: pre_step_action => rk23_pre_step
            !! Performs any pre-step actions.
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

    type, extends(single_step_integrator) :: runge_kutta_853
        !! An 8th order Dormand-Prince type 8th order integrator.
        real(real64), private, allocatable, dimension(:) :: k2
        real(real64), private, allocatable, dimension(:) :: k3
        real(real64), private, allocatable, dimension(:) :: k4
        real(real64), private, allocatable, dimension(:) :: k5
        real(real64), private, allocatable, dimension(:) :: k6
        real(real64), private, allocatable, dimension(:) :: k7
        real(real64), private, allocatable, dimension(:) :: k8
        real(real64), private, allocatable, dimension(:) :: k9
        real(real64), private, allocatable, dimension(:) :: k10
        real(real64), private, allocatable, dimension(:) :: yerr2
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
        real(real64), private, allocatable, dimension(:) :: rc4
        real(real64), private, allocatable, dimension(:) :: rc5
        real(real64), private, allocatable, dimension(:) :: rc6
        real(real64), private, allocatable, dimension(:) :: rc7
        real(real64), private, allocatable, dimension(:) :: rc8
        real(real64), private, allocatable, dimension(:) :: work
        real(real64), private :: m_stepSize
    contains
        procedure, public :: pre_step_action => rk853_pre_step
            !! Performs any pre-step actions.
        procedure, private :: initialize => rk853_init
            !! Initializes private storage arrays for the integrator.
        procedure, private :: initialize_interp => rk853_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: get_order => rk853_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk853_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: attempt_step => rk853_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rk853_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rk853_interp
            !! Performs the interpolation.
        procedure, public :: compute_error_norm => rk853_estimate_error
            !! Computes the norm of the scaled error estimate.
    end type

    type, abstract, extends(single_step_integrator) :: &
        diagonally_implicit_runge_kutta
        !! Defines a diagonally implicit Runge-Kutta (DIRK) type integrator.
        !!
        !! Remarks:
        !!
        !! The integrators based upon this type are expected to utilize a
        !! constant on their diagonal such that \( a_{ii} = \gamma \) thereby
        !! allowing for a constant iteration matrix.  If the method does not
        !! employ such behaviors, it is recommended to overload the Newton 
        !! solve routine and code the appropriate behavior.
        integer(int32), private :: m_maxNewtonIter = 20
            ! Maximum allowable number of Newton iterations.
        real(real64), private :: m_newtonTol = 1.0d-6
            ! Newton iteration convergence tolerance.
    contains
        procedure, public :: get_max_newton_iteration_count => &
            dirk_get_max_newton_iter
            !! Gets the maximum number of Newton iterations allowed.
        procedure, public :: set_max_newton_iteration_count => &
            dirk_set_max_newton_iter
            !! Sets the maximum number of Newton iterations allowed.
        procedure, public :: get_newton_tolerance => dirk_get_newton_tol
            !! Gets the convergence tolerance for the Newton iteration.
        procedure, public :: set_newton_tolerance => dirk_set_newton_tol
            !! Sets the convergence tolerance for the Newton iteration.
    end type
contains
! ******************************************************************************
! RUNGE_KUTTA_45
! ------------------------------------------------------------------------------
subroutine rk45_pre_step(this, prevs, sys, h, x, y, f, err)
    !! Placeholder routine for any pre-step actions.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    logical, intent(in) :: prevs
        !! Defines the status of the previous step.  The value is true 
        !! if the previous step was successful; else, false if the 
        !! previous step failed.
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Process
    return
end subroutine

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

    yn = y + h * a21 * f
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
subroutine rk45_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn)
    use diffeq_dprk45_constants
    !! Sets up the interpolation process.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
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

    ! Quick Return
    if (.not.dense) return

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
subroutine rk23_pre_step(this, prevs, sys, h, x, y, f, err)
    !! Placeholder routine for any pre-step actions.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    logical, intent(in) :: prevs
        !! Defines the status of the previous step.  The value is true 
        !! if the previous step was successful; else, false if the 
        !! previous step failed.
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Process
    return
end subroutine

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

    yn = y + h * a21 * f
    call sys%fcn(x + h * c2, yn, this%k2)

    yn = y + h * (a31 * f + a32 * this%k2)
    call sys%fcn(x + h * c3, yn, this%k3)

    yn = y + h * (a41 * f + a42 * this%k2 + a43 * this%k3)
    call sys%fcn(x + h * c4, yn, fn)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * this%k2 + e3 * this%k3 + e4 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk32_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn)
    use diffeq_bsrk32_constants
    !! Sets up the interpolation process.
    class(runge_kutta_23), intent(inout) :: this
        !! The runge_kutta_23 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
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

    ! Quick Return
    if (.not.dense) return

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

! ******************************************************************************
! RUNGE_KUTTA_853
! ------------------------------------------------------------------------------
subroutine rk853_pre_step(this, prevs, sys, h, x, y, f, err)
    !! Placeholder routine for any pre-step actions.
    class(runge_kutta_853), intent(inout) :: this
        !! The runge_kutta_853 object.
    logical, intent(in) :: prevs
        !! Defines the status of the previous step.  The value is true 
        !! if the previous step was successful; else, false if the 
        !! previous step failed.
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Process
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine rk853_init(this, neqn)
    !! Initializes private storage arrays for the integrator.
    class(runge_kutta_853), intent(inout) :: this
        !! The runge_kutta_853 object.
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
            deallocate(this%k7)
            deallocate(this%k8)
            deallocate(this%k9)
            deallocate(this%k10)
            deallocate(this%yerr2)
        end if
    end if

    ! If we're here, we need to allocate the storage arrays
    allocate( &
        this%k2(neqn), &
        this%k3(neqn), &
        this%k4(neqn), &
        this%k5(neqn), &
        this%k6(neqn), &
        this%k7(neqn), &
        this%k8(neqn), &
        this%k9(neqn), &
        this%k10(neqn), &
        this%yerr2(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk853_init_interp(this, neqn)
    !! Allocates storage for the interpolation process.
    class(runge_kutta_853), intent(inout) :: this
        !! The runge_kutta_853 object.
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
            deallocate(this%rc6)
            deallocate(this%rc7)
            deallocate(this%rc8)
            deallocate(this%work)
        end if
    end if

    allocate( &
        this%rc1(neqn), &
        this%rc2(neqn), &
        this%rc3(neqn), &
        this%rc4(neqn), &
        this%rc5(neqn), &
        this%rc6(neqn), &
        this%rc7(neqn), &
        this%rc8(neqn), &
        this%work(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
pure function rk853_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(runge_kutta_853), intent(in) :: this
        !! The runge_kutta_853 object.
    integer(int32) :: rst
        !! The order.
    rst = 8
end function

! ------------------------------------------------------------------------------
pure function rk853_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(runge_kutta_853), intent(in) :: this
        !! The runge_kutta_853 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine rk853_attempt_step(this, sys, h, x, y, f, yn, fn, yerr)
    use diffeq_rk853_constants
    !! Attempts an integration step for this integrator.
    class(runge_kutta_853), intent(inout) :: this
        !! The runge_kutta_853 object.
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
    ! k1 = f as the derivatives were computed outside of this routine

    yn = y + h * a21
    call sys%fcn(x + c2 * h, yn, this%k2)

    yn = y + h * (a31 * f + a32 * this%k2)
    call sys%fcn(x + c3 * h, yn, this%k3)

    yn = y + h * (a41 * f + a43 * this%k3)
    call sys%fcn(x + c4 * h, yn, this%k4)

    yn = y + h * (a51 * f + a53 * this%k3 + a54 * this%k4)
    call sys%fcn(x + c5 * h, yn, this%k5)

    yn = y + h * (a61 * f + a64 * this%k4 + a65 * this%k5)
    call sys%fcn(x + c6 * h, yn, this%k6)

    yn = y + h * (a71 * f + a74 * this%k4 + a75 * this%k5 + a76 * this%k6)
    call sys%fcn(x + c7 * h, yn, this%k7)

    yn = y + h * (a81 * f + a84 * this%k4 + a85 * this%k5 + a86 * this%k6 + &
        a87 * this%k7)
    call sys%fcn(x + c8 * h, yn, this%k8)

    yn = y + h * (a91 * f + a94 * this%k4 + a95 * this%k5 + a96 * this%k6 + &
        a97 * this%k7 + a98 * this%k8)
    call sys%fcn(x + c9 * h, yn, this%k9)

    yn = y + h * (a101 * f + a104 * this%k4 + a105 * this%k5 + &
        a106 * this%k6 + a107 * this%k7 + a108 * this%k8 + a109 * this%k9)
    call sys%fcn(x + c10 * h, yn, this%k10)

    yn = y + h * (a111 * f + a114 * this%k4 + a115 * this%k5 + &
        a116 * this%k6 + a117 * this%k7 + a118 * this%k8 + a119 * this%k9 + &
        a1110 * this%k10)
    call sys%fcn(x + c11 * h, yn, this%k2)

    yn = y + h * (a121 * f + a124 * this%k4 + a125 * this%k5 + &
        a126 * this%k6 + a127 * this%k7 + a128 * this%k8 + a129 * this%k9 + &
        a1210 * this%k10 + a1211 * this%k2)
    call sys%fcn(x + h, yn, this%k3)

    this%k4 = b1 * f + b6 * this%k6 + b7 * this%k7 + b8 * this%k8 + &
        b9 * this%k9 + b10 * this%k10 + b11 * this%k2 + b12 * this%k3
    yn = y + h * this%k4
    call sys%fcn(x + h, yn, fn)

    yerr = this%k4 - bhh1 * f - bhh2 * this%k9 - bhh3 * this%k3
    this%yerr2 = er1 * f + er6 * this%k6 + er7 * this%k7 + er8 * this%k8 + &
        er9 * this%k9 + er10 * this%k10 + er11 * this%k2 + er12 * this%k3
    this%m_stepSize = h
end subroutine

! ------------------------------------------------------------------------------
subroutine rk853_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn)
    use diffeq_rk853_constants
    !! Sets up the interpolation process.
    class(runge_kutta_853), intent(inout) :: this
        !! The runge_kutta_853 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
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
    real(real64) :: ydiff, bspl, h

    ! Quick Return
    if (.not.dense) return

    ! Initialization
    n = size(y)
    h = xn - x
    call this%initialize_interp(n)

    ! Process
    do i = 1, n
        this%rc1(i) = y(i)
        ydiff = yn(i) - y(i)
        this%rc2(i) = ydiff
        bspl = h * f(i) - ydiff
        this%rc3(i) = bspl
        this%rc4(i) = ydiff - h * fn(i) - bspl
        this%rc5(i) = d41 * f(i) + d46 * this%k6(i) + d47 * this%k7(i) + &
            d48 * this%k8(i) + d49 * this%k9(i) + d410 * this%k10(i) + &
            d411 * this%k2(i) + d412 * this%k3(i)
        this%rc6(i) = d51 * f(i) + d56 * this%k6(i) + d57 * this%k7(i) + &
            d58 * this%k8(i) + d59 * this%k9(i) + d510 * this%k10(i) + &
            d511 * this%k2(i) + d512 * this%k3(i)
        this%rc7(i) = d61 * f(i) + d66 * this%k6(i) + d67 * this%k7(i) + &
            d68 * this%k8(i) + d69 * this%k9(i) + d610 * this%k10(i) + &
            d611 * this%k2(i) + d612 * this%k3(i)
        this%rc8(i) = d71 * f(i) + d76 * this%k6(i) + d77 * this%k7(i) + &
            d78 * this%k8(i) + d79 * this%k9(i) + d710 * this%k10(i) + &
            d711 * this%k2(i) + d712 * this%k3(i)
    end do

    this%work = y + h * (a141 * f + a147 * this%k7 + a148 * this%k8 + &
        a149 * this%k9 + a1410 * this%k10 + a1411 * this%k2 + &
        a1412 * this%k3 + a1413 * fn)
    call sys%fcn(x + c14 * h, this%work, this%k10)

    this%work = y + h * (a151 * f + a156 * this%k6 + a157 * this%k7 + &
        a158 * this%k8 + a1511 * this%k2 + a1512 * this%k3 + a1513 * fn + &
        a1514 * this%k10)
    call sys%fcn(x + c15 * h, this%work, this%k2)

    this%work = y + h * (a161 * f + a166 * this%k6 + a167 * this%k7 + &
        a168 * this%k8 + a169 * this%k9 + a1613 * fn + a1614 * this%k10 + &
        a1615 * this%k2)
    call sys%fcn(x + c16 * h, this%work, this%k3)

    do i = 1, n
        this%rc5(i) = h * (this%rc5(i) + d413 * fn(i) + d414 * this%k10(i) + &
            d415 * this%k2(i) + d416 * this%k3(i))
        this%rc6(i) = h * (this%rc6(i) + d513 * fn(i) + d514 * this%k10(i) + &
            d515 * this%k2(i) + d516 * this%k3(i))
        this%rc7(i) = h * (this%rc7(i) + d613 * fn(i) + d614 * this%k10(i) + &
            d615 * this%k2(i) + d616 * this%k3(i))
        this%rc8(i) = h * (this%rc8(i) + d713 * fn(i) + d714 * this%k10(i) + &
            d715 * this%k2(i) + d716 * this%k3(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine rk853_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(runge_kutta_853), intent(in) :: this
        !! The runge_kutta_853 object.
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
    y = this%rc1 + s * (this%rc2 + s1 * (this%rc3 + s * (this%rc4 + s1 * &
        (this%rc5 + s * (this%rc6 + s1 * (this%rc7 + s * this%rc8))))))
end subroutine

! ------------------------------------------------------------------------------
pure function rk853_estimate_error(this, y, yest, yerr) result(rst)
    !! Computes the norm of the scaled error estimate.  A value less than one
    !! indicates a successful step.  A value greater than one suggests that the
    !! results do not meet the requested tolerances.
    class(runge_kutta_853), intent(in) :: this
        !! The runge_kutta_853 object.
    real(real64), intent(in), dimension(:) :: y
        !! The previously accepted solution array (N-element).
    real(real64), intent(in), dimension(size(y)) :: yest
        !! An N-element array containing the next solution point estimate.
    real(real64), intent(in), dimension(size(y)) :: yerr
        !! An N-element array containing the estimate of error for each
        !! equation.
    real(real64) :: rst
        !! The norm of the scaled error.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: err, err2, sf, denom, atol, rtol

    ! Initialization
    n = size(y)
    atol = this%get_absolute_tolerance()
    rtol = this%get_relative_tolerance()
    err = 0.0d0
    err2 = 0.0d0
    
    ! Process
    do i = 1, n
        sf = atol + rtol * max(abs(y(i)), abs(yest(i)))
        err2 = err2 + (yerr(i) / sf)**2
        err = err + (this%yerr2(i) / sf)**2
    end do
    denom = err + 1.0d-2 * err2
    if (denom <= 0.0d0) denom = 1.0d0
    rst = err * sqrt(1.0d0 / (n * denom)) * abs(this%m_stepSize)
end function

! ******************************************************************************
! DIAGONALLY IMPLICIT INTEGRATOR
! ------------------------------------------------------------------------------
pure function dirk_get_max_newton_iter(this) result(rst)
    !! Gets the maximum number of Newton iterations allowed.
    class(diagonally_implicit_runge_kutta), intent(in) :: this
        !! The diagonally_implicit_runge_kutta object.
    integer(int32) :: rst
        !! The iteration limit.
    rst = this%m_maxNewtonIter
end function

! --------------------
subroutine dirk_set_max_newton_iter(this, x)
    !! Sets the maximum number of Newton iterations allowed.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
    integer(int32), intent(in) :: x
        !! The iteration limit
    this%m_maxNewtonIter = x
end subroutine

! ------------------------------------------------------------------------------
pure function dirk_get_newton_tol(this) result(rst)
    !! Gets the convergence tolerance for the Newton iteration.
    class(diagonally_implicit_runge_kutta), intent(in) :: this
        !! The diagonally_implicit_runge_kutta object.
    real(real64) :: rst
        !! The tolerance.
    rst = this%m_newtonTol
end function

! --------------------
subroutine dirk_set_newton_tol(this, x)
    !! Sets the convergence tolerance for the Newton iteration.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
    real(real64), intent(in) :: x
        !! The tolerance.
    this%m_newtonTol = x
end subroutine

! ------------------------------------------------------------------------------
subroutine dirk_solve_newton(this, i, sys, h, x, y)
    !! Solves the Newton iteration problem for the i-th stage assuming a
    !! the coefficient matrix is constant on its diagonal such that 
    !! \( a_{ii} = \gamma \).  This constraint allows for a constant iteration
    !! matrix.
    class(diagonally_implicit_runge_kutta), intent(in) :: this
        !! The diagonally_implicit_runge_kutta object.
    integer(int32), intent(in) :: i
        !! The number of the stage to solve.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the dependent
        !! variables.
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module