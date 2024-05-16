module diffeq_runge_kutta
    use iso_fortran_env
    use diffeq_errors
    use diffeq_base
    use ferror
    implicit none
    private
    public :: runge_kutta_45
    public :: runge_kutta_23
    public :: runge_kutta_853
    public :: diagonally_implicit_runge_kutta
    public :: dirk_array_value
    public :: dirk_matrix_value
    public :: dirk_integer_inquiry
    public :: dirk_action
    public :: implicit_runge_kutta_4

    type, extends(single_step_integrator) :: runge_kutta_45
        !! The Dormand-Prince, Runge-Kutta integrator (4th/5th order).
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
        procedure, public :: get_stage_count => rk45_get_stage_count
            !! Gets the stage count for this integrator.
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
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
    contains
        procedure, public :: pre_step_action => rk23_pre_step
            !! Performs any pre-step actions.
        procedure, private :: initialize_interp => rk23_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: get_order => rk23_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk23_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => rk23_get_stage_count
            !! Gets the stage count for this integrator.
        procedure, public :: attempt_step => rk23_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rk23_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rk23_interp
            !! Performs the interpolation.
    end type

    type, extends(single_step_integrator) :: runge_kutta_853
        !! An 8th order Dormand-Prince type 8th order integrator.
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
        ! procedure, private :: initialize => rk853_init
            !! Initializes private storage arrays for the integrator.
        procedure, private :: initialize_interp => rk853_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: get_order => rk853_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk853_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => rk853_get_stage_count
            !! Gets the stage count for this integrator.
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
        !! employ such behaviors, it is recommended to overload the  
        !! newton_iteration routine and code the appropriate behavior.
        integer(int32), private :: m_maxNewtonIter = 20
            ! Maximum allowable number of Newton iterations.
        real(real64), private :: m_newtonTol = 1.0d-6
            ! Newton iteration convergence tolerance.
        real(real64), private, allocatable, dimension(:,:) :: a
            ! The Newton iteration matrix.
        integer(int32), private, allocatable, dimension(:) :: pvt
            ! The LU-factorization pivot tracking array.
        real(real64), private, allocatable, dimension(:,:) :: jac
            ! The Jacobian matrix.
        real(real64), private, allocatable, dimension(:,:) :: mass
            ! The mass matrix.
        logical, private :: m_updateJacobian = .true.
            ! True if the Jacobian should be updated; else, false.
    contains
        procedure(dirk_matrix_value), public, pass, deferred :: &
            get_model_coefficient
            !! Gets the requested model coefficient from the Butcher tableau.
        procedure(dirk_array_value), public, pass, deferred :: get_weight
            !! Gets a weighting value from the Butcher tableau.
        procedure(dirk_array_value), public, pass, deferred :: get_node
            !! Gets a node value from the Butcher tableau.
        procedure(dirk_array_value), public, pass, deferred :: &
            get_error_coefficient
            !! Gets the requested coefficient of the error stage.
        procedure(dirk_integer_inquiry), public, pass, deferred :: &
            get_stage_count
            !! Gets the stage count for the integrator.
        procedure(dirk_action), public, pass, deferred :: fill_table
            !! Fills out the Butcher tableau for this integrator object.
        procedure, private :: initialize_matrices => dirk_init_matrices
            !! Allocates internal storage for the system matrices.
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
        procedure, public :: newton_iteration => dirk_solve_newton
            !! Solves the Newton iteration problem for the i-th stage.
        procedure, public :: pre_step_action => dirk_form_matrix
            !! Constructs the system matrix.
        procedure, public :: attempt_step => dirk_attempt_step
            !! Attempts an integration step for this integrator.
    end type

    interface
        pure function dirk_matrix_value(this, i, j) result(rst)
            !! Gets a value from a stored matrix.
            use iso_fortran_env
            import diagonally_implicit_runge_kutta
            class(diagonally_implicit_runge_kutta), intent(in) :: this
                !! The diagonally_implicit_runge_kutta object.
            integer(int32), intent(in) :: i
                !! The row index.
            integer(int32), intent(in) :: j
                !! The column index.
            real(real64) :: rst
                !! The requested value.
        end function

        pure function dirk_array_value(this, i) result(rst)
            !! Gets a value from a stored array.
            use iso_fortran_env
            import diagonally_implicit_runge_kutta
            class(diagonally_implicit_runge_kutta), intent(in) :: this
                !! The diagonally_implicit_runge_kutta object.
            integer(int32), intent(in) :: i
                !! The index.
            real(real64) :: rst
                !! The requested value.
        end function

        pure function dirk_integer_inquiry(this) result(rst)
            !! Gets an integer value from the diagonally_implicit_runge_kutta
            !! object.
            use iso_fortran_env
            import diagonally_implicit_runge_kutta
            class(diagonally_implicit_runge_kutta), intent(in) :: this
                !! The diagonally_implicit_runge_kutta object.
            integer(int32) :: rst
                !! The requested value.
        end function

        subroutine dirk_action(this)
            !! Performs an action on the integrator object.
            import diagonally_implicit_runge_kutta
            class(diagonally_implicit_runge_kutta), intent(inout) :: this
                !! The diagonally_implicit_runge_kutta object.
        end subroutine
    end interface

    type, extends(diagonally_implicit_runge_kutta) :: implicit_runge_kutta_4
        !! A singly-diagonally implicit 4th order Runge-Kutta integrator
        !! suitable for integrating stiff systems of differential equations.
        real(real64), private :: m_coeffs(6,6)
            ! Butcher tableau coefficients.
        real(real64), private :: m_weights(6)
            ! Butcher tableau weighting factors.
        real(real64), private :: m_nodes(6)
            ! Butcher tableau nodes.
        real(real64), private :: m_errorCoeffs(6)
            ! Error estimator coefficients.
        real(real64), private :: m_interpCoeffs(4,6)
            ! Interpolation constants
        logical, private :: m_filledTable = .false.
            ! True if the tableau variables have been populated; else, false.
        real(real64), private, allocatable, dimension(:,:) :: m_k
            ! An N-by-NSTAGES matrix containing the derivatives at each stage.
    contains
        procedure, public :: fill_table => irk4_fill_coeffs
            !! Fills out the Butcher tableau for this integrator object.
        procedure, public :: get_order => irk4_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => irk4_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => irk4_get_stage_count
            !! Gets the number of stages for this integrator.
        procedure, public :: post_step_action => irk4_set_up_interp
            !! Sets up the interpolation process.
        procedure, public :: interpolate => irk4_interp
            !! Performs the interpolation.
        procedure, public :: get_model_coefficient => irk4_get_model_coeff
            !! Gets the requested model coefficient from the Butcher tableau.
        procedure, public :: get_weight => irk4_get_weight
            !! Gets a weighting value from the Butcher tableau.
        procedure, public :: get_node => irk4_get_node
            !! Gets a node value from the Butcher tableau.
        procedure, public :: get_error_coefficient => irk4_get_error_coeff
            !! Gets the requested coefficient of the error stage.
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
pure function rk45_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    integer(int32) :: rst
        !! The stage count.
    rst = 7
end function

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
subroutine rk45_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k)
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
    real(real64), intent(out), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

    ! Local Variables
    integer(int32) :: n

    ! Initialization
    n = size(y)

    ! Process
    ! k1 = f as the derivatives were computed from the previous step
    k(:,1) = f

    yn = y + h * a21 * f
    call sys%fcn(x + h * c2, yn, k(:,2))

    yn = y + h * (a31 * f + a32 * k(:,2))
    call sys%fcn(x + h * c3, yn, k(:,3))

    yn = y + h * (a41 * f + a42 * k(:,2) + a43 * k(:,3))
    call sys%fcn(x + h * c4, yn, k(:,4))

    yn = y + h * (a51 * f + a52 * k(:,2) + a53 * k(:,3) + a54 * k(:,4))
    call sys%fcn(x + h * c5, yn, k(:,5))

    yn = y + h * (a61 * f + a62 * k(:,2) + a63 * k(:,3) + a64 * k(:,4) + &
        a65 * k(:,5))
    call sys%fcn(x + h * c6, yn, k(:,6))

    yn = y + h * (a71 * f + a73 * k(:,3) + a74 * k(:,4) + a75 * k(:,5) + a76 * k(:,6))
    call sys%fcn(x + h * c7, yn, fn)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * k(:,2) + e3 * k(:,3) + e4 * k(:,4) + &
        e5 * k(:,5) + e6 * k(:,6) + e7 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k)
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
    real(real64), intent(inout), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

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
        this%rc5 = h * (d1 * f(i) + d3 * k(i,3) + d4 * k(i,4) + &
            d5 * k(i,5) + d6 * k(i,6) + d7 * fn(i))
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
pure function rk23_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    integer(int32) :: rst
        !! The order.
    rst = 3
end function

! ------------------------------------------------------------------------------
pure function rk23_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function rk23_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(runge_kutta_23), intent(in) :: this
        !! The runge_kutta_23 object.
    integer(int32) :: rst
        !! The stage count.
    rst = 4
end function

! ------------------------------------------------------------------------------
subroutine rk23_init_interp(this, neqn)
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
subroutine rk23_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k)
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
    real(real64), intent(out), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

    ! Local Variables
    integer(int32) :: n

    ! Initialization
    n = size(y)

    ! Process
    ! k1 = f as the derivatives were computed from the previous step
    k(:,1) = f

    yn = y + h * a21 * f
    call sys%fcn(x + h * c2, yn, k(:,2))

    yn = y + h * (a31 * f + a32 * k(:,2))
    call sys%fcn(x + h * c3, yn, k(:,3))

    yn = y + h * (a41 * f + a42 * k(:,2) + a43 * k(:,3))
    call sys%fcn(x + h * c4, yn, fn)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * k(:,2) + e3 * k(:,3) + e4 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk23_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k)
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
    real(real64), intent(inout), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

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
subroutine rk23_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
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
pure function rk853_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(runge_kutta_853), intent(in) :: this
        !! The runge_kutta_853 object.
    integer(int32) :: rst
        !! The stage count.
    rst = 12
end function

! ------------------------------------------------------------------------------
subroutine rk853_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k)
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
    real(real64), intent(out), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

    ! Local Variables
    integer(int32) :: n

    ! Initialization
    n = size(y)
    if (.not.allocated(this%yerr2)) allocate(this%yerr2(n))

    ! Process
    ! k1 = f as the derivatives were computed outside of this routine
    k(:,1) = f

    yn = y + h * a21
    call sys%fcn(x + c2 * h, yn, k(:,2))

    yn = y + h * (a31 * f + a32 * k(:,2))
    call sys%fcn(x + c3 * h, yn, k(:,3))

    yn = y + h * (a41 * f + a43 * k(:,3))
    call sys%fcn(x + c4 * h, yn, k(:,4))

    yn = y + h * (a51 * f + a53 * k(:,3) + a54 * k(:,4))
    call sys%fcn(x + c5 * h, yn, k(:,5))

    yn = y + h * (a61 * f + a64 * k(:,4) + a65 * k(:,5))
    call sys%fcn(x + c6 * h, yn, k(:,6))

    yn = y + h * (a71 * f + a74 * k(:,4) + a75 * k(:,5) + a76 * k(:,6))
    call sys%fcn(x + c7 * h, yn, k(:,7))

    yn = y + h * (a81 * f + a84 * k(:,4) + a85 * k(:,5) + a86 * k(:,6) + &
        a87 * k(:,7))
    call sys%fcn(x + c8 * h, yn, k(:,8))

    yn = y + h * (a91 * f + a94 * k(:,4) + a95 * k(:,5) + a96 * k(:,6) + &
        a97 * k(:,7) + a98 * k(:,8))
    call sys%fcn(x + c9 * h, yn, k(:,9))

    yn = y + h * (a101 * f + a104 * k(:,4) + a105 * k(:,5) + &
        a106 * k(:,6) + a107 * k(:,7) + a108 * k(:,8) + a109 * k(:,9))
    call sys%fcn(x + c10 * h, yn, k(:,10))

    yn = y + h * (a111 * f + a114 * k(:,4) + a115 * k(:,5) + &
        a116 * k(:,6) + a117 * k(:,7) + a118 * k(:,8) + a119 * k(:,9) + &
        a1110 * k(:,10))
    call sys%fcn(x + c11 * h, yn, k(:,2))

    yn = y + h * (a121 * f + a124 * k(:,4) + a125 * k(:,5) + &
        a126 * k(:,6) + a127 * k(:,7) + a128 * k(:,8) + a129 * k(:,9) + &
        a1210 * k(:,10) + a1211 * k(:,2))
    call sys%fcn(x + h, yn, k(:,3))

    k(:,4) = b1 * f + b6 * k(:,6) + b7 * k(:,7) + b8 * k(:,8) + &
        b9 * k(:,9) + b10 * k(:,10) + b11 * k(:,2) + b12 * k(:,3)
    yn = y + h * k(:,4)
    call sys%fcn(x + h, yn, fn)

    yerr = k(:,4) - bhh1 * f - bhh2 * k(:,9) - bhh3 * k(:,3)
    this%yerr2 = er1 * f + er6 * k(:,6) + er7 * k(:,7) + er8 * k(:,8) + &
        er9 * k(:,9) + er10 * k(:,10) + er11 * k(:,2) + er12 * k(:,3)
    this%m_stepSize = h
end subroutine

! ------------------------------------------------------------------------------
subroutine rk853_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k)
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
    real(real64), intent(inout), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

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
        this%rc5(i) = d41 * f(i) + d46 * k(i,6) + d47 * k(i,7) + &
            d48 * k(i,8) + d49 * k(i,9) + d410 * k(i,10) + &
            d411 * k(i,2) + d412 * k(i,3)
        this%rc6(i) = d51 * f(i) + d56 * k(i,6) + d57 * k(i,7) + &
            d58 * k(i,8) + d59 * k(i,9) + d510 * k(i,10) + &
            d511 * k(i,2) + d512 * k(i,3)
        this%rc7(i) = d61 * f(i) + d66 * k(i,6) + d67 * k(i,7) + &
            d68 * k(i,8) + d69 * k(i,9) + d610 * k(i,10) + &
            d611 * k(i,2) + d612 * k(i,3)
        this%rc8(i) = d71 * f(i) + d76 * k(i,6) + d77 * k(i,7) + &
            d78 * k(i,8) + d79 * k(i,9) + d710 * k(i,10) + &
            d711 * k(i,2) + d712 * k(i,3)
    end do

    this%work = y + h * (a141 * f + a147 * k(:,7) + a148 * k(:,8) + &
        a149 * k(:,9) + a1410 * k(:,10) + a1411 * k(:,2) + &
        a1412 * k(:,3) + a1413 * fn)
    call sys%fcn(x + c14 * h, this%work, k(:,10))

    this%work = y + h * (a151 * f + a156 * k(:,6) + a157 * k(:,7) + &
        a158 * k(:,8) + a1511 * k(:,2) + a1512 * k(:,3) + a1513 * fn + &
        a1514 * k(:,10))
    call sys%fcn(x + c15 * h, this%work, k(:,2))

    this%work = y + h * (a161 * f + a166 * k(:,6) + a167 * k(:,7) + &
        a168 * k(:,8) + a169 * k(:,9) + a1613 * fn + a1614 * k(:,10) + &
        a1615 * k(:,2))
    call sys%fcn(x + c16 * h, this%work, k(:,3))

    do i = 1, n
        this%rc5(i) = h * (this%rc5(i) + d413 * fn(i) + d414 * k(i,10) + &
            d415 * k(i,2) + d416 * k(i,3))
        this%rc6(i) = h * (this%rc6(i) + d513 * fn(i) + d514 * k(i,10) + &
            d515 * k(i,2) + d516 * k(i,3))
        this%rc7(i) = h * (this%rc7(i) + d613 * fn(i) + d614 * k(i,10) + &
            d615 * k(i,2) + d616 * k(i,3))
        this%rc8(i) = h * (this%rc8(i) + d713 * fn(i) + d714 * k(i,10) + &
            d715 * k(i,2) + d716 * k(i,3))
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
subroutine dirk_solve_newton(this, i, sys, h, x, y, f, niter, &
    success, err)
    use linalg, only : solve_lu
    !! Solves the Newton iteration problem for the i-th stage assuming a
    !! the coefficient matrix is constant on its diagonal such that 
    !! \( a_{ii} = \gamma \).  This constraint allows for a constant iteration
    !! matrix.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
    integer(int32), intent(in) :: i
        !! The number of the stage to solve.
    class(ode_container), intent(inout) :: sys
        !! The N-element pivot tracking array from the LU factorization
        !! process.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(inout), dimension(:) :: y
        !! An N-element array containing the initial guess at the solution.
        !! On output, the updated solution estimate.
    real(real64), intent(inout), dimension(:) :: f
        !! An N-element array, that on input, contains the values of the 
        !! derivatives at x.  On output, the updated derivative estimate.
    integer(int32), intent(out) :: niter
        !! The number of iterations performed.
    logical, intent(out) :: success
        !! Returns true if the iteration process successfully converged; else,
        !! false if the iteration could not converge.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
        !!      memory allocation issue.

    ! Local Variables
    integer(int32) :: j, n, flag, maxiter
    real(real64) :: tol, z, alpha
    real(real64), allocatable, dimension(:) :: w, dy
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(y)
    maxiter = this%get_max_newton_iteration_count()
    tol = this%get_newton_tolerance()
    alpha = this%get_model_coefficient(i, i)
    success = .false.
    allocate(w(n), dy(n), stat = flag, source = 0.0d0)
    if (flag /= 0) then
        call report_memory_error(errmgr, "dirk_solve_newton", flag)
        return
    end if

    ! Ensure the tableau is filled
    call this%fill_table()

    ! Process
    do j = 1, i - 1
        w = w + this%get_model_coefficient(i,j) * f
    end do
    w = y + h * w
    z = x + this%get_node(i) * h
    call sys%fcn(z, w, f)

    ! Iteration
    do niter = 1, maxiter
        ! Compute the right-hand-side
        dy = w + h * alpha * f - y

        ! Solve the system
        call solve_lu(this%a, this%pvt, dy)

        ! Update the solution
        y = y + dy
        call sys%fcn(z, y, f)

        ! Check for convergence
        if (norm2(dy) < tol) then
            success = .true.
            exit
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine dirk_init_matrices(this, n, usemass)
    !! Allocates internal storage for the system matrices.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
    integer(int32), intent(in) :: n
        !! The number of equations being integrated.
    logical, intent(in) :: usemass
        !! True if a mass matrix is used; else, false.

    ! Local Variables
    integer(int32) :: ns

    ! Process
    ns = this%get_stage_count()
    if (allocated(this%jac)) then
        if (size(this%jac, 1) == n .and. size(this%jac, 2) == n) then
            ! All is good
            return
        else
            deallocate(this%jac)
            if (allocated(this%mass)) deallocate(this%mass)
            deallocate(this%pvt)
            deallocate(this%a)
        end if
    end if
    if (usemass) then
        allocate( &
            this%jac(n, n), &
            this%mass(n, n), &
            this%pvt(n), &
            this%a(n, n) &
        )
    else
        allocate( &
            this%jac(n, n), &
            this%pvt(n), &
            this%a(n, n) &
        )
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine dirk_form_matrix(this, prevs, sys, h, x, y, f, err)
    use linalg, only : lu_factor
    !! Constructs the system matrix of the form \( A = f M - J \), and then 
    !! computes it's LU factorization.  The LU-factored form of A is stored 
    !! internally.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
    logical, intent(in) :: prevs
        !! Defines the status of the previous step.  The value is true if the
        !! previous step was successful; else, false if the previous step 
        !! failed.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives at x.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Local Variables
    integer(int32) :: i, n, ns
    logical :: useMass
    real(real64) :: fac
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Quick Return
    if (.not.this%m_updateJacobian) return

    ! Ensure the tableau is filled
    call this%fill_table()
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(y)
    useMass = associated(sys%mass_matrix)
    call this%initialize_matrices(n, useMass)
    ns = this%get_stage_count()
    fac = this%get_model_coefficient(ns, ns) * h

    ! Process
    call sys%compute_jacobian(x, y, this%jac, errmgr)
    if (errmgr%has_error_occurred()) return
    if (useMass) then
        this%a = this%mass - fac * this%jac
    else
        this%a = -fac * this%jac
        do i = 1, n
            this%a(i,i) = this%a(i,i) + 1.0d0
        end do
    end if

    ! Compute the LU factorization
    call lu_factor(this%a, this%pvt)

    ! The Jacobian is updated
    this%m_updateJacobian = .false.
end subroutine
! ------------------------------------------------------------------------------
subroutine dirk_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k)
    !! Attempts an integration step for this integrator.
    class(diagonally_implicit_runge_kutta), intent(inout) :: this
        !! The diagonally_implicit_runge_kutta object.
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
    real(real64), intent(out), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

    ! Local Variables
    logical :: success
    integer(int32) :: i, n, ns, niter, maxiter, itertracking

    ! Initialization
    n = size(y)
    ns = this%get_stage_count()
    maxiter = this%get_max_newton_iteration_count()

    ! Ensure the tableau is filled
    call this%fill_table()

    ! Cycle over each stage and solve the Newton problem
    itertracking = 0
    yn = y  ! use the previously accepted solution as an initial guess
    k(:,1) = f ! along with the corresponding derivatives
    do i = 1, ns
        call this%newton_iteration(i, sys, h, x, yn, k(:,i), niter, &
            success)
        itertracking = max(niter, itertracking)
        if (.not.success) exit
    end do

    ! Do we need to update the Jacobian?
    if (.not.success) then
        ! No convergence could be achieved - force a Jacobian update
        this%m_updateJacobian = .true.

        ! Do we need to update the step size?
    end if

    ! Update the solution estimate
    do i = 1, ns
        if (i == 1) then
            yn = this%get_weight(i) * k(:,i)
            yerr = this%get_error_coefficient(i) * k(:,i)
        else
            yn = yn + this%get_weight(i) * k(:,i)
            yerr = yerr + this%get_error_coefficient(i) * k(:,i)
        end if
    end do
    yn = y + h * yn
    yerr = h * yerr

    ! Do we need to update the Jacobian
    if (itertracking > maxiter / 2) then
        ! We're having a tough enough time to justify a change
        this%m_updateJacobian = .true.
    end if
end subroutine

! ******************************************************************************
! 4th ORDER IMPLICIT RUNGE KUTTA
! ------------------------------------------------------------------------------
pure function irk4_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32) :: rst
        !! The order.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function irk4_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function


! ------------------------------------------------------------------------------
pure function irk4_get_stage_count(this) result(rst)
    !! Gets the stage count for the integrator.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32) :: rst
        !! The stage count.
    rst = 6
end function

! ------------------------------------------------------------------------------
subroutine irk4_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k)
    !! Sets up the interpolation process.
    class(implicit_runge_kutta_4), intent(inout) :: this
        !! The implicit_runge_kutta_4 object.
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
    real(real64), intent(inout), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.

    ! Local Variables
    integer(int32) :: n, nstages

    ! Initialization
    n = size(y)
    nstages = this%get_stage_count()
    if (allocated(this%m_k)) then
        ! Ensure it's sized appropriately
        if (size(this%m_k, 1) /= n .or. size(this%m_k, 2) /= nstages) then
            deallocate(this%m_k)
            allocate(this%m_k(n, nstages))
        end if
    else
        allocate(this%m_k(n, nstages))
    end if
    this%m_k = k    ! store the derivative values for use by the interpolation
end subroutine

! ------------------------------------------------------------------------------
subroutine irk4_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
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
    integer(int32) :: i, j, norder, nstages, n
    real(real64) :: h, theta, bi
    real(real64), allocatable, dimension(:) :: w

    ! Initialization
    n = size(yn)
    norder = this%get_order()
    nstages = this%get_stage_count()
    h = xn1 - xn
    theta = (x - xn) / h
    allocate(w(n), source = 0.0d0)

    ! Process
    do i = 1, nstages
        bi = 0.0d0
        do j = 1, norder
            bi = bi + this%m_interpCoeffs(j,i) * theta**j
        end do
        w = w + bi * this%m_k(:,i)
    end do
    y = yn + h * w
end subroutine

! ------------------------------------------------------------------------------
subroutine irk4_fill_coeffs(this)
    use diffeq_sdirk4_constants
    !! Fills the internal arrays storing the model coefficients.
    class(implicit_runge_kutta_4), intent(inout) :: this
        !! The implicit_runge_kutta_4 object.

    ! Quick Return
    if (this%m_filledTable) return

    ! A
    this%m_coeffs = 0.0d0

    this%m_coeffs(2,1) = a21
    this%m_coeffs(2,2) = a22

    this%m_coeffs(3,1) = a31
    this%m_coeffs(3,2) = a32
    this%m_coeffs(3,3) = a33

    this%m_coeffs(4,1) = a41
    this%m_coeffs(4,2) = a42
    this%m_coeffs(4,3) = a43
    this%m_coeffs(4,4) = a44

    this%m_coeffs(5,1) = a51
    this%m_coeffs(5,2) = a52
    this%m_coeffs(5,3) = a53
    this%m_coeffs(5,4) = a54
    this%m_coeffs(5,5) = a55

    this%m_coeffs(6,1) = b1
    this%m_coeffs(6,2) = b2
    this%m_coeffs(6,3) = b3
    this%m_coeffs(6,4) = b4
    this%m_coeffs(6,5) = b5
    this%m_coeffs(6,6) = b6

    ! B
    this%m_weights(1) = b1
    this%m_weights(2) = b2
    this%m_weights(3) = b3
    this%m_weights(4) = b4
    this%m_weights(5) = b5
    this%m_weights(6) = b6

    ! C
    this%m_nodes(1) = 0.0d0
    this%m_nodes(2) = c2
    this%m_nodes(3) = c3
    this%m_nodes(4) = c4
    this%m_nodes(5) = c5
    this%m_nodes(6) = c6

    ! E
    this%m_errorCoeffs(1) = b1a - b1
    this%m_errorCoeffs(2) = b2a - b2
    this%m_errorCoeffs(3) = b3a - b3
    this%m_errorCoeffs(4) = b4a - b4
    this%m_errorCoeffs(5) = b5a - b5
    this%m_errorCoeffs(6) = b6a - b6

    ! Interpolation Constants
    this%m_interpCoeffs(1,1) = bs11
    this%m_interpCoeffs(2,1) = bs21
    this%m_interpCoeffs(3,1) = bs31
    this%m_interpCoeffs(4,1) = bs41

    this%m_interpCoeffs(1,2) = bs12
    this%m_interpCoeffs(2,2) = bs22
    this%m_interpCoeffs(3,2) = bs32
    this%m_interpCoeffs(4,2) = bs42

    this%m_interpCoeffs(1,3) = bs13
    this%m_interpCoeffs(2,3) = bs23
    this%m_interpCoeffs(3,3) = bs33
    this%m_interpCoeffs(4,3) = bs43

    this%m_interpCoeffs(1,4) = bs14
    this%m_interpCoeffs(2,4) = bs24
    this%m_interpCoeffs(3,4) = bs34
    this%m_interpCoeffs(4,4) = bs44

    this%m_interpCoeffs(1,5) = bs15
    this%m_interpCoeffs(2,5) = bs25
    this%m_interpCoeffs(3,5) = bs35
    this%m_interpCoeffs(4,5) = bs45

    this%m_interpCoeffs(1,6) = bs16
    this%m_interpCoeffs(2,6) = bs26
    this%m_interpCoeffs(3,6) = bs36
    this%m_interpCoeffs(4,6) = bs46

    ! Update status
    this%m_filledTable = .true.
end subroutine

! ------------------------------------------------------------------------------
pure function irk4_get_model_coeff(this, i, j) result(rst)
    !! Gets the requested model coefficient from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    integer(int32), intent(in) :: j
        !! The column index.
    real(real64) :: rst
        !! The coefficient

    ! Return the requested value
    rst = this%m_coeffs(i, j)
end function

! ------------------------------------------------------------------------------
pure function irk4_get_weight(this, i) result(rst)
    !! Gets the requested weighting factor from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    real(real64) :: rst
        !! The weighting factor.

    ! Return the requested value
    rst = this%m_weights(i)
end function

! ------------------------------------------------------------------------------
pure function irk4_get_node(this, i) result(rst)
    !! Gets the requested node from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    real(real64) :: rst
        !! The node value.

    ! Return the requested value
    rst = this%m_nodes(i)
end function

! ------------------------------------------------------------------------------
pure function irk4_get_error_coeff(this, i) result(rst)
    !! Gets the requested error estimator coefficient from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    real(real64) :: rst
        !! The coefficient.

    ! Return the requested value
    rst = this%m_errorCoeffs(i)
end function

! ------------------------------------------------------------------------------
end module