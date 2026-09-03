module diffeq_multistep
    !! Provides the base structure for implicit, variable step-size, variable
    !! order, multi-step integrators.
    !!
    !! A multi-step method advances a system written as
    !! \[
    !! M(x,y) y' = f(x,y)
    !! \]
    !! by relating the solution at \(x_{n+1}\) to several previously accepted
    !! solution points rather than to intermediate stage values.  The routines
    !! collected here supply the machinery common to such methods: storage and
    !! management of the solution history, Newton iteration on the implicit
    !! corrector equation, polynomial interpolation of the history for both
    !! prediction and dense output, and the adaptive step-size and order
    !! controller.
    !!
    !! A concrete integrator supplies only the method-specific pieces by
    !! implementing predictor, correct, get_maximum_supported_order,
    !! get_error_constant, and get_truncation_constant.  Methods whose
    !! corrector equation is solved by Newton iteration should extend
    !! newton_multistep_integrator, which supplies that corrector.
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use linalg, only : solve_least_squares_full
    implicit none
    private
    public :: multistep_integrator
    public :: newton_multistep_integrator
    public :: build_linearized_system
    public :: predictor_step
    public :: multistep_corrector
    public :: multistep_integer_inquiry
    public :: multistep_error_constant
    public :: newton_divided_differences
    public :: evaluate_newton_polynomial

    real(real64), parameter :: min_error = 1.0d-10
        ! A floor applied to scaled error norms so that step size ratios remain
        ! finite when a step is essentially exact.

    type, extends(ode_integrator), abstract :: multistep_integrator
        !! Defines the structure of a multi-step integrator.
        integer(int32), private :: m_maxOrder = 0
            ! The user-requested upper limit on the order of the method.  A
            ! value below one deprives the request of meaning and is taken to
            ! mean that the limit of the method itself should govern.
        integer(int32), private :: m_order = 1
            ! The order currently in use.
        integer(int32), private :: m_stepsAtOrder = 0
            ! The number of consecutive successful steps taken at m_order.
        integer(int32), private :: m_nHist = 0
            ! The number of solution points currently stored in the history.
        real(real64), private, allocatable, dimension(:) :: m_xHist
            ! Independent variable history, most recent value first.
        real(real64), private, allocatable, dimension(:,:) :: m_yHist
            ! Solution history, most recent value in column one.
        real(real64), private, allocatable, dimension(:,:) :: m_fHist
            ! Derivative history, most recent value in column one.
    contains
        procedure(predictor_step), public, pass, deferred :: predictor
            !! Evaluates the predictor.
        procedure(multistep_corrector), public, pass, deferred :: correct
            !! Refines the prediction by applying the implicit formula.
        procedure(multistep_integer_inquiry), public, pass, deferred :: &
            get_maximum_supported_order
            !! Gets the highest order the method is permitted to use.
        procedure(multistep_error_constant), public, pass, deferred :: &
            get_error_constant
            !! Gets the error constant for the requested order.
        procedure(multistep_error_constant), public, pass, deferred :: &
            get_truncation_constant
            !! Gets the truncation error constant for the requested order.
        procedure, public :: evaluate_derivative => mi_evaluate_derivative
            !! Evaluates the derivative stored in the solution history.
        procedure, public :: prepare_for_solve => mi_prepare_for_solve
            !! Prepares internal state for the beginning of an integration.
        procedure, public :: get_order => mi_get_order
            !! Gets the order currently in use.
        procedure, public :: get_maximum_order => mi_get_max_order
            !! Gets the upper limit on the order of the method.
        procedure, public :: set_maximum_order => mi_set_max_order
            !! Sets the upper limit on the order of the method.
        procedure, public :: get_stored_step_count => mi_get_stored_count
            !! Gets the number of solution points held in the history.
        procedure, public :: get_stored_variable => mi_get_stored_x
            !! Gets a stored independent variable value.
        procedure, public :: get_stored_solution => mi_get_stored_y
            !! Gets a stored solution vector.
        procedure, public :: get_stored_derivative => mi_get_stored_f
            !! Gets a stored derivative vector.
        procedure, public :: reset_history => mi_reset_history
            !! Clears the stored solution history.
        procedure, public :: interpolate => mi_interpolate
            !! Interpolates the solution history to the requested point.
        procedure, public :: solve => mi_ode_solver
            !! Solves the supplied system of ODE's.
        procedure, private :: initialize_history => mi_init_history
            !! Allocates storage for the solution history.
        procedure, private :: push_history => mi_push_history
            !! Stores a newly accepted solution point.
        procedure, private :: effective_max_order => mi_effective_max_order
            !! Gets the order limit honoring both user and method limits.
        procedure, private :: adjust_order_and_step => mi_adjust_order_step
            !! Selects the order and step size for the next step.
    end type

    type, extends(multistep_integrator), abstract :: newton_multistep_integrator
        !! Defines the structure of a multi-step integrator whose corrector
        !! equation is solved by Newton iteration.
        !!
        !! Solving the corrector by Newton iteration requires a Jacobian, but
        !! in return it places no restriction upon the stiffness of the problem
        !! and permits a singular mass matrix, as the mass matrix then appears
        !! only as a factor within the iteration matrix.
        integer(int32), private :: m_newtonIterLimit = 12
            ! The limit on the number of Newton iterations allowed.
        real(real64), private :: m_newtonTolerance = 1.0d-8
            ! The convergence tolerance on the Euclidean norm of the residual
            ! for the Newton iteration process.
    contains
        procedure(build_linearized_system), public, pass, deferred :: &
            build_linear_system
            !! Builds the linearized system for the implicit corrector.
        procedure, public :: correct => mi_newton_iteration
            !! Performs the Newton iteration process.
        procedure, public :: newton_iteration => mi_newton_iteration
            !! Performs the Newton iteration process.
        procedure, public :: get_newton_step_limit => mi_get_newton_step_limit
            !! Gets the maximum number of Newton iterations allowed.
        procedure, public :: set_newton_step_limit => mi_set_newton_step_limit
            !! Sets the maximum number of Newton iterations allowed.
        procedure, public :: get_newton_tolerance => mi_get_newton_tolerance
            !! Gets the Newton iteration convergence tolerance.
        procedure, public :: set_newton_tolerance => mi_set_newton_tolerance
            !! Sets the Newton iteration convergence tolerance.
    end type

    interface
        subroutine build_linearized_system(this, sys, h, x, y, A, b, args)
            !! Builds the linearized system used to compute the solution for the
            !! implicit step via Newton iteration.  This routine forms matrices
            !! \(A\) and \(b\) of the system \(A \delta = b\).
            use iso_fortran_env, only : real64
            import newton_multistep_integrator
            import ode_container
            class(newton_multistep_integrator), intent(inout) :: this
                !! The newton_multistep_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! The current value(s) of the dependent variable(s).
            real(real64), intent(out), dimension(:,:) :: A
                !! The N-by-N matrix \(A\), where, for BDF methods at least,
                !! \(A = \frac{\alpha_o}{h_n} M - J\).
            real(real64), intent(out), dimension(:) :: b
                !! The N-element right-hand-side vector for the linear system.
                !! This vector is the negative of the corrector residual so
                !! that the Newton update is \(y \leftarrow y + \delta\).
            class(*), intent(inout), optional, target :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        subroutine predictor_step(this, sys, h, x, y, f, yn, fn, args)
            !! Evaluates the predictor.
            !!
            !! Besides supplying the initial estimate for the Newton iteration,
            !! this routine is the natural place for a concrete integrator to
            !! compute any step-dependent quantities, such as method
            !! coefficients, that the corrector requires.  It is called exactly
            !! once per step attempt.
            use iso_fortran_env, only : real64
            import multistep_integrator
            import ode_container
            class(multistep_integrator), intent(inout) :: this
                !! The multistep_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! The current value(s) of the dependent variable(s).
            real(real64), intent(in), dimension(:) :: f
                !! An N-element array containing the values of the derivatives
                !! at x.  For a system carrying a mass matrix these are the
                !! values of \(y'\) rather than of \(f\).
            real(real64), intent(out), dimension(:) :: yn
                !! An N-element array where this routine will write the next
                !! solution estimate at x + h.
            real(real64), intent(out), dimension(:) :: fn
                !! An N-element array where this routine will write the next
                !! derivative estimate at x + h.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        subroutine multistep_corrector(this, sys, h, x, ypred, ynew, &
            converged, args)
            !! Refines the prediction by applying the implicit formula.
            use iso_fortran_env, only : real64
            import multistep_integrator
            import ode_container
            class(multistep_integrator), intent(inout) :: this
                !! The multistep_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The current step size.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: ypred
                !! The N-element array containing the predicted values at x + h.
            real(real64), intent(out), dimension(:) :: ynew
                !! The new estimate of the dependent variables.
            logical, intent(out) :: converged
                !! True if the corrector succeeded; else, false, in which case
                !! the step is retried on a shorter step size.
            class(*), intent(inout), optional, target :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        pure function multistep_integer_inquiry(this) result(rst)
            !! Returns an integer value from the multistep_integrator object.
            use iso_fortran_env, only : int32
            import multistep_integrator
            class(multistep_integrator), intent(in) :: this
                !! The multistep_integrator object.
            integer(int32) :: rst
                !! The requested value.
        end function

        pure function multistep_error_constant(this, k) result(rst)
            !! Returns the constant relating the difference between the
            !! corrected and predicted solutions to the local truncation error
            !! of the method at the requested order.
            use iso_fortran_env, only : int32, real64
            import multistep_integrator
            class(multistep_integrator), intent(in) :: this
                !! The multistep_integrator object.
            integer(int32), intent(in) :: k
                !! The order of interest.
            real(real64) :: rst
                !! The error constant.
        end function
    end interface
contains
! ------------------------------------------------------------------------------
pure function mi_get_newton_step_limit(this) result(rst)
    !! Gets the maximum number of Newton iterations allowed.
    class(newton_multistep_integrator), intent(in) :: this
        !! The newton_multistep_integrator object.
    integer(int32) :: rst
        !! The iteration limit.
    rst = this%m_newtonIterLimit
end function

! --------------------
subroutine mi_set_newton_step_limit(this, x)
    !! Sets the maximum number of Newton iterations allowed.
    class(newton_multistep_integrator), intent(inout) :: this
        !! The newton_multistep_integrator object.
    integer(int32), intent(in) :: x
        !! The iteration limit.
    this%m_newtonIterLimit = x
end subroutine

! ------------------------------------------------------------------------------
pure function mi_get_newton_tolerance(this) result(rst)
    !! Gets the convergence tolerance on the Euclidean norm of the residual
    !! for the Newton iteration process.
    class(newton_multistep_integrator), intent(in) :: this
        !! The newton_multistep_integrator object.
    real(real64) :: rst
        !! The convergence tolerance.
    rst = this%m_newtonTolerance
end function

! --------------------
subroutine mi_set_newton_tolerance(this, x)
    !! Sets the convergence tolerance on the Euclidean norm of the residual
    !! for the Newton iteration process.
    class(newton_multistep_integrator), intent(inout) :: this
        !! The newton_multistep_integrator object.
    real(real64), intent(in) :: x
        !! The convergence tolerance.
    this%m_newtonTolerance = x
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_newton_iteration(this, sys, h, x, ypred, ynew, converged, args)
    !! Performs the Newton iteration process.
    !!
    !! The iteration solves the implicit corrector equation supplied by
    !! build_linear_system.  That system is solved by means of a
    !! rank-revealing orthogonal factorization, so a singular mass matrix is
    !! tolerated provided the iteration matrix itself remains nonsingular, as
    !! is the case for an index-1 problem.
    class(newton_multistep_integrator), intent(inout) :: this
        !! The newton_multistep_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: ypred
        !! The N-element array containing the predicted values at x + h.
    real(real64), intent(out), dimension(:) :: ynew
        !! The new estimate of the dependent variables.
    logical, intent(out) :: converged
        !! True if the Newton iteration process converged; else, false.
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information in and 
        !! out of the differential equation subroutine.

    ! Local Variables
    integer(int32) :: i, n, nmax
    real(real64) :: tol
    real(real64), allocatable, dimension(:,:) :: A
    real(real64), allocatable, dimension(:) :: b, delta

    ! Initialization
    converged = .false.
    n = size(ypred)
    allocate(A(n,n), b(n))
    nmax = this%get_newton_step_limit()
    tol = this%get_newton_tolerance()
    ynew = ypred

    ! Input Check
    if (size(ynew) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR

    ! Process
    do i = 1, nmax
        ! Create the system matrices
        call this%build_linear_system(sys, h, x + h, ynew, A, b, args)

        ! Solve the system
        delta = solve_least_squares_full(A, b)

        ! Update the solution and test for convergence.  The tolerance is
        ! applied relative to the size of the solution so that the test
        ! behaves consistently across problems of differing scale.
        ynew = ynew + delta
        if (norm2(delta) <= tol * max(1.0d0, norm2(ynew))) then
            converged = .true.
            exit
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_evaluate_derivative(this, sys, x, y, dydx, args)
    !! Evaluates the derivative that the solution history stores.
    !!
    !! The default implementation simply evaluates the right-hand side, which
    !! is correct whenever the mass matrix is the identity.  An integrator that
    !! requires \(y'\) itself, rather than \(f\), overrides this routine to
    !! resolve the mass matrix.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x
        !! The value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution at x.
    real(real64), intent(out), dimension(:) :: dydx
        !! An N-element array where the derivative will be written.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information in and out
        !! of the differential equation subroutine.
    call sys%fcn(x, y, dydx, args)
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_prepare_for_solve(this)
    !! Performs whatever preparation the integrator requires at the beginning
    !! of an integration.
    !!
    !! The default implementation does nothing.  An integrator that caches a
    !! matrix factorization between steps overrides this routine to discard it,
    !! so that reusing the object with a different model cannot pick up stale
    !! information.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
end subroutine

! ------------------------------------------------------------------------------
pure function mi_get_order(this) result(rst)
    !! Gets the order of the integrator currently in use.
    !!
    !! A variable-order method changes this value as the integration proceeds,
    !! so it reports the order of the most recently attempted step.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32) :: rst
        !! The order.
    rst = this%m_order
end function

! ------------------------------------------------------------------------------
pure function mi_get_max_order(this) result(rst)
    !! Gets the upper limit on the order of the method.
    !!
    !! The value reported is the limit actually in force, which is the lesser
    !! of any limit requested by the caller and the highest order the method
    !! supports.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32) :: rst
        !! The order limit.
    rst = this%effective_max_order()
end function

! --------------------
subroutine mi_set_max_order(this, x)
    !! Sets the upper limit on the order of the method.
    !!
    !! The value is constrained to the range supported by the method.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    integer(int32), intent(in) :: x
        !! The order limit.
    if (x < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_maxOrder = min(x, this%get_maximum_supported_order())
end subroutine

! --------------------
pure function mi_effective_max_order(this) result(rst)
    !! Gets the order limit honoring both the user-requested limit and the
    !! highest order the method itself supports.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32) :: rst
        !! The order limit.
    if (this%m_maxOrder < 1) then
        rst = this%get_maximum_supported_order()
    else
        rst = min(this%m_maxOrder, this%get_maximum_supported_order())
    end if
    rst = max(1, rst)
end function

! ------------------------------------------------------------------------------
pure function mi_get_stored_count(this) result(rst)
    !! Gets the number of previously accepted solution points currently held in
    !! the history.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32) :: rst
        !! The number of stored points.
    rst = this%m_nHist
end function

! --------------------
pure function mi_get_stored_x(this, i) result(rst)
    !! Gets a stored value of the independent variable, where an index of one
    !! denotes the most recently accepted point.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32), intent(in) :: i
        !! The one-based index into the history.
    real(real64) :: rst
        !! The stored value.
    rst = this%m_xHist(i)
end function

! --------------------
pure function mi_get_stored_y(this, i) result(rst)
    !! Gets a stored solution vector, where an index of one denotes the most
    !! recently accepted point.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32), intent(in) :: i
        !! The one-based index into the history.
    real(real64), allocatable, dimension(:) :: rst
        !! The stored solution.
    rst = this%m_yHist(:,i)
end function

! --------------------
pure function mi_get_stored_f(this, i) result(rst)
    !! Gets a stored derivative vector, where an index of one denotes the most
    !! recently accepted point.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    integer(int32), intent(in) :: i
        !! The one-based index into the history.
    real(real64), allocatable, dimension(:) :: rst
        !! The stored derivative.
    rst = this%m_fHist(:,i)
end function

! ------------------------------------------------------------------------------
subroutine mi_reset_history(this)
    !! Clears the stored solution history.
    !!
    !! The integrator restarts at first order the next time it is used.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    this%m_nHist = 0
    this%m_order = 1
    this%m_stepsAtOrder = 0
end subroutine

! --------------------
subroutine mi_init_history(this, neqn)
    !! Allocates storage for the solution history.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.

    ! Local Variables
    integer(int32) :: nmax

    ! Room is required for the highest usable order plus the additional points
    ! the order controller needs in order to test one order higher.
    nmax = this%effective_max_order() + 3
    if (allocated(this%m_xHist)) then
        if (size(this%m_yHist, 1) == neqn .and. &
            size(this%m_xHist) == nmax) return
        deallocate(this%m_xHist)
        deallocate(this%m_yHist)
        deallocate(this%m_fHist)
    end if
    allocate(this%m_xHist(nmax), source = 0.0d0)
    allocate(this%m_yHist(neqn, nmax), source = 0.0d0)
    allocate(this%m_fHist(neqn, nmax), source = 0.0d0)
end subroutine

! --------------------
subroutine mi_push_history(this, x, y, f)
    !! Stores a newly accepted solution point at the front of the history.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    real(real64), intent(in) :: x
        !! The independent variable value.
    real(real64), intent(in), dimension(:) :: y
        !! The solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! The derivative at x.

    ! Local Variables
    integer(int32) :: i, nmax

    ! Process
    nmax = size(this%m_xHist)
    do i = min(this%m_nHist, nmax - 1), 1, -1
        this%m_xHist(i+1) = this%m_xHist(i)
        this%m_yHist(:,i+1) = this%m_yHist(:,i)
        this%m_fHist(:,i+1) = this%m_fHist(:,i)
    end do
    this%m_xHist(1) = x
    this%m_yHist(:,1) = y
    this%m_fHist(:,1) = f
    this%m_nHist = min(this%m_nHist + 1, nmax)
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_interpolate(this, x, y)
    !! Interpolates the stored solution history to the requested value of the
    !! independent variable.
    !!
    !! The polynomial passing through the most recently accepted points is the
    !! same polynomial on which the method itself is based, so dense output
    !! carries the accuracy of the current order and costs no additional
    !! function evaluations.
    class(multistep_integrator), intent(in) :: this
        !! The multistep_integrator object.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which to compute the
        !! interpolation.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where the interpolated solution is written.

    ! Local Variables
    integer(int32) :: np
    real(real64), allocatable, dimension(:,:) :: dd

    ! Input Checking
    if (this%m_nHist < 1) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (size(y) /= size(this%m_yHist, 1)) error stop DIFFEQ_ARRAY_SIZE_ERROR

    ! Process
    np = min(this%m_order + 1, this%m_nHist)
    if (np == 1) then
        y = this%m_yHist(:,1)
        return
    end if
    allocate(dd(size(y), np))
    call newton_divided_differences(this%m_xHist(1:np), &
        this%m_yHist(:,1:np), dd)
    call evaluate_newton_polynomial(this%m_xHist(1:np), dd, x, y)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine newton_divided_differences(x, y, d)
    !! Computes the Newton divided-difference coefficients of the polynomial
    !! interpolating the supplied points.
    !!
    !! The resulting coefficients define
    !! \[
    !! P(t) = d_1 + d_2 (t - x_1) + d_3 (t - x_1)(t - x_2) + \cdots
    !! \]
    !! Divided differences are symmetric in their arguments, so the nodes need
    !! not be sorted.  Storing them most-recent-first keeps the terms of the
    !! extrapolation in decreasing order of significance.
    real(real64), intent(in), dimension(:) :: x
        !! An M-element array containing the interpolation nodes.
    real(real64), intent(in), dimension(:,:) :: y
        !! An N-by-M matrix whose columns contain the values at each node.
    real(real64), intent(out), dimension(:,:) :: d
        !! An N-by-M matrix where the coefficients will be written.

    ! Local Variables
    integer(int32) :: i, j, m

    ! Process
    m = size(x)
    d(:,1:m) = y(:,1:m)
    do j = 1, m - 1
        do i = m, j + 1, -1
            d(:,i) = (d(:,i) - d(:,i-1)) / (x(i) - x(i-j))
        end do
    end do
end subroutine

! --------------------
pure subroutine evaluate_newton_polynomial(x, d, xi, y)
    !! Evaluates a polynomial stored in Newton divided-difference form.
    real(real64), intent(in), dimension(:) :: x
        !! An M-element array containing the interpolation nodes.
    real(real64), intent(in), dimension(:,:) :: d
        !! An N-by-M matrix containing the divided-difference coefficients.
    real(real64), intent(in) :: xi
        !! The value at which to evaluate the polynomial.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where the result will be written.

    ! Local Variables
    integer(int32) :: i, m

    ! Process
    m = size(x)
    y = d(:,m)
    do i = m - 1, 1, -1
        y = d(:,i) + (xi - x(i)) * y
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_adjust_order_step(this, yprev, ynew, ek, h)
    !! Selects the order and step size to use for the next step.
    !!
    !! The order is chosen by asking which of the neighboring orders would
    !! permit the largest step, and the step size then follows from the error
    !! estimate belonging to the order selected.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    real(real64), intent(in), dimension(:) :: yprev
        !! The solution at the beginning of the step just accepted.
    real(real64), intent(in), dimension(:) :: ynew
        !! The solution at the end of the step just accepted.
    real(real64), intent(in) :: ek
        !! The scaled error norm of the step just accepted.
    real(real64), intent(inout) :: h
        !! On input, the step size just used.  On output, the step size to
        !! attempt next.

    ! Local Variables
    integer(int32) :: k, kmax, kbest, mmax, n
    real(real64) :: ekc, ekm1, ekp1, r, rbest, safe
    real(real64), allocatable, dimension(:) :: errv
    real(real64), allocatable, dimension(:,:) :: dd

    ! Initialization
    n = size(ynew)
    k = this%m_order
    kmax = this%effective_max_order()
    safe = this%get_step_size_factor()
    ekc = ek
    ekm1 = -1.0d0
    ekp1 = -1.0d0

    ! Build the divided-difference table over the accepted history.  The table
    ! supplies the local error the method would have committed had a
    ! neighboring order been used, and that is what drives the order selection.
    ! Every candidate is measured the same way so that the comparison between
    ! them is not distorted by the differing constants of each order.
    mmax = min(this%m_nHist, k + 3)
    if (mmax >= 2) then
        allocate(dd(n, mmax))
        call newton_divided_differences(this%m_xHist(1:mmax), &
            this%m_yHist(:,1:mmax), dd)

        if (mmax >= k + 2) then
            errv = this%get_truncation_constant(k) * &
                history_difference(this%m_xHist, dd, k + 1)
            ekc = this%compute_error_norm(yprev, ynew, errv)
        end if

        if (k > 1 .and. mmax >= k + 1) then
            errv = this%get_truncation_constant(k - 1) * &
                history_difference(this%m_xHist, dd, k)
            ekm1 = this%compute_error_norm(yprev, ynew, errv)
        end if

        if (k < kmax .and. mmax >= k + 3) then
            errv = this%get_truncation_constant(k + 1) * &
                history_difference(this%m_xHist, dd, k + 2)
            ekp1 = this%compute_error_norm(yprev, ynew, errv)
        end if
    end if

    ! Select whichever order permits the largest step
    kbest = k
    rbest = step_ratio(ekc, k)
    if (ekm1 >= 0.0d0) then
        r = step_ratio(ekm1, k - 1)
        if (r > rbest) then
            rbest = r
            kbest = k - 1
        end if
    end if

    ! An increase in order is considered only once the current order has proven
    ! itself, which keeps the method from chasing transients.
    if (ekp1 >= 0.0d0 .and. this%m_stepsAtOrder >= k + 1) then
        r = step_ratio(ekp1, k + 1)
        if (r > rbest) then
            rbest = r
            kbest = k + 1
        end if
    end if

    ! A formula of order k relates k + 1 points, and the predictor that
    ! accompanies it must extrapolate an equal number, so the order is held to
    ! one less than the number of points on hand.
    kbest = max(1, min(kbest, kmax, this%m_nHist - 1))
    if (kbest == this%m_order) then
        this%m_stepsAtOrder = this%m_stepsAtOrder + 1
    else
        this%m_order = kbest
        this%m_stepsAtOrder = 0
    end if

    ! Apply the step size ratio.  Modest changes are suppressed because the
    ! stability of a multi-step formula degrades when the spacing of the
    ! history it interpolates is disturbed without cause.
    r = safe * rbest
    if (r < 1.0d0) then
        r = max(r, 1.0d-1)
    else if (r < 2.0d0) then
        r = 1.0d0
    else
        r = min(r, 1.0d1)
    end if
    h = r * h

contains
    pure function step_ratio(e, korder) result(rst)
        !! Estimates the factor by which the step size could change were the
        !! given order used.
        real(real64), intent(in) :: e
        integer(int32), intent(in) :: korder
        real(real64) :: rst
        rst = (1.0d0 / max(e, min_error))**(1.0d0 / (korder + 1.0d0))
    end function

    pure function history_difference(xh, dtab, m) result(rst)
        !! Forms the difference between the accepted solution and the value a
        !! degree m-1 extrapolation of the preceding points would have given.
        real(real64), intent(in), dimension(:) :: xh
        real(real64), intent(in), dimension(:,:) :: dtab
        integer(int32), intent(in) :: m
        real(real64), allocatable, dimension(:) :: rst
        integer(int32) :: i
        real(real64) :: p
        p = 1.0d0
        do i = 2, m + 1
            p = p * (xh(1) - xh(i))
        end do
        rst = dtab(:,m+1) * p
    end function
end subroutine

! ------------------------------------------------------------------------------
subroutine mi_ode_solver(this, sys, x, iv, args)
    !! Solves the supplied system of ODE's.
    !!
    !! The integration is driven by a predictor-corrector cycle.  Each step
    !! extrapolates the solution history to form a prediction, refines that
    !! prediction with Newton iteration on the implicit formula, and compares
    !! the two to estimate the local error.  A step whose error satisfies the
    !! requested tolerances is accepted and appended to the history; otherwise
    !! it is retried on a shorter step.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array of independent variable values at which to return the
        !! solution to the ODE's.  If only two values are supplied, the
        !! solution is returned at every successful integration step.
    real(real64), intent(in), dimension(:) :: iv
        !! An array containing the initial values for each ODE.
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    logical :: dense, converged
    integer(int32) :: i, j, n, neqn, nsteps, k, nfail
    real(real64) :: h, xo, xn, xmax, e, minstep
    real(real64), allocatable, dimension(:) :: f, y, yn, fn, fnew, ypred, &
        yerr, yi

    ! Initialization
    n = size(x)

    ! Input Checking
    if (n < 2) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (.not.sys%get_is_ode_defined()) error stop DIFFEQ_MISSING_ARGUMENT_ERROR

    ! Additional Initialization
    neqn = size(iv)
    xo = x(1)
    xmax = x(n)
    j = 2
    dense = (n > 2)
    nsteps = this%get_step_limit()
    minstep = this%get_minimum_step_size()
    nfail = 0

    ! Memory Allocations
    allocate( &
        f(neqn), &
        y(neqn), &
        yn(neqn), &
        fn(neqn), &
        fnew(neqn), &
        ypred(neqn), &
        yerr(neqn), &
        yi(neqn) &
    )

    ! The method restarts at first order.  A first-order formula needs no
    ! history beyond the initial condition, so no starting procedure and no
    ! consistent initial derivative are required.
    call this%reset_history()
    call this%prepare_for_solve()
    call this%initialize_history(neqn)
    y = iv

    ! Estimate an initial step size
    !
    ! Outputs:
    ! - f: Value of the derivatives at xo
    ! - h: Initial step size estimate
    call this%estimate_inital_step_size(sys, xo, xmax, iv, f, h, args)

    ! Seed the history
    call this%evaluate_derivative(sys, xo, y, f, args)
    call this%push_history(xo, y, f)

    ! Store the initial conditions
    call this%append_to_buffer(xo, y)

    ! Cycle until integration is complete
    do i = 1, nsteps
        ! The order can never exceed what the stored history will support
        this%m_order = max(1, min(this%m_order, this%m_nHist - 1))
        k = this%m_order

        ! Do not overshoot the end of the integration range if the solver has
        ! been asked to terminate on it exactly
        if (.not.this%get_allow_overshoot()) then
            if (abs(xo + h) > abs(xmax)) h = xmax - xo
        end if

        ! Predict, then correct
        call this%predictor(sys, h, xo, y, f, ypred, fn, args)
        call this%correct(sys, h, xo, ypred, yn, converged, args)

        if (.not.converged) then
            ! A failure of the corrector is a statement about the step size
            ! rather than about the accuracy of the result, so retry on a
            ! shorter step and fall back to the most robust order if failures
            ! persist.
            nfail = nfail + 1
            this%m_stepsAtOrder = 0
            h = 2.5d-1 * h
            if (nfail >= 3) this%m_order = 1
            if (abs(h) < minstep) error stop DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
            cycle
        end if

        ! The difference between the corrected and predicted solutions is
        ! proportional to the local truncation error of the formula
        yerr = this%get_error_constant(k) * (yn - ypred)
        e = this%compute_error_norm(y, yn, yerr)

        ! The negated test also rejects a step whose error is not a number,
        ! which is how a non-stiff corrector announces that it has diverged.
        if (.not.(e <= 1.0d0)) then
            ! Reject the step and try again on a shorter one
            nfail = nfail + 1
            this%m_stepsAtOrder = 0
            if (e <= huge(e)) then
                h = h * max(1.0d-1, min(5.0d-1, this%get_step_size_factor() * &
                    (1.0d0 / max(e, min_error))**(1.0d0 / (k + 1.0d0))))
            else
                h = 2.5d-1 * h
            end if
            if (nfail >= 3) this%m_order = 1
            if (abs(h) < minstep) error stop DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
            cycle
        end if

        ! The step was successful
        nfail = 0
        xn = xo + h
        call this%evaluate_derivative(sys, xn, yn, fnew, args)
        call this%push_history(xn, yn, fnew)

        ! Store the results
        if (dense) then
            ! Interpolate to each requested point falling within this step
            interp : do while (abs(x(j)) <= abs(xn))
                call this%interpolate(x(j), yi)
                call this%append_to_buffer(x(j), yi)
                j = j + 1
                if (j > n) exit interp
            end do interp
        else if (abs(xn) > abs(xmax)) then
            ! The step carried past the end of the range, so report the
            ! solution interpolated back onto the requested end point
            call this%interpolate(xmax, yi)
            call this%append_to_buffer(xmax, yi)
        else
            ! Report every successful step
            call this%append_to_buffer(xn, yn)
        end if

        ! Are we done?
        if (abs(xn) >= abs(xmax)) return

        ! Choose the order and step size for the next step
        call this%adjust_order_and_step(y, yn, e, h)
        if (abs(h) > this%get_maximum_step_size()) then
            h = sign(this%get_maximum_step_size(), h)
        end if
        if (abs(h) < minstep) error stop DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR

        ! Update for the next step
        xo = xn
        y = yn
        f = fnew
    end do

    ! If we're here, the solver has run out of allowable steps
    error stop DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR
end subroutine

! ------------------------------------------------------------------------------
end module