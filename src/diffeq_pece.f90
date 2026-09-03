module diffeq_pece
    !! Explicit predictor, implicit corrector (PECE) multi-step solvers.
    !!
    !! An Adams method advances the solution by integrating an interpolating
    !! polynomial fitted to the derivative history,
    !! \[
    !! y_{n+1} = y_n + \int_{x_n}^{x_{n+1}} P(t) \, dt.
    !! \]
    !! The predictor fits that polynomial to derivatives already known, giving
    !! an explicit Adams-Bashforth formula, while the corrector includes the
    !! derivative at the new point as well, giving an implicit Adams-Moulton
    !! formula of the same order.  Evaluating the two in sequence, and taking
    !! their difference as a measure of the local error, is the classical PECE
    !! cycle.  Each step therefore costs exactly two derivative evaluations and
    !! requires neither a Jacobian nor the solution of a linear system.
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use diffeq_multistep
    use lapack, only : DGETRF, DGETRS
    implicit none
    private
    public :: adams

    integer(int32), parameter :: adams_order_limit = 12
        ! The order is capped at twelve so that the tabulated error constants
        ! below span every order the method is permitted to select.

    type, extends(multistep_integrator) :: adams
        !! A variable step-size, variable order (1-12) Adams-Bashforth-Moulton
        !! integrator operating in predictor-corrector (PECE) mode.
        !!
        !! The corrector is applied once rather than iterated to convergence,
        !! which is what distinguishes PECE from a fully implicit method.  A
        !! single correction is only contractive while \(h L < 1\), where
        !! \(L\) bounds the Lipschitz constant of the problem, so this
        !! integrator is intended for non-stiff systems.  Applied to a stiff
        !! problem it holds its accuracy only by driving the step size down to
        !! what stability permits, which can cost orders of magnitude more
        !! steps than an implicit method requires; such a problem belongs to
        !! bdf or to one of the implicit Runge-Kutta integrators.
        !!
        !! A system carrying a mass matrix is supported, but the matrix must be
        !! nonsingular.  The method needs \(y'\) itself, so the mass matrix is
        !! factored and resolved rather than merely appearing as a factor, and
        !! a singular matrix is reported as an error.  A differential-algebraic
        !! problem must therefore be given to bdf.
        real(real64), private, allocatable, dimension(:) :: m_fpred
            ! The derivative evaluated at the predicted solution.
        real(real64), private, allocatable, dimension(:,:) :: m_lu
            ! The LU factorization of the mass matrix.
        integer(int32), private, allocatable, dimension(:) :: m_pivot
            ! The pivot array of the mass matrix factorization.
        logical, private :: m_massFactored = .false.
            ! True if the mass matrix factorization is current.
    contains
        procedure, public :: predictor => adams_predictor
            !! Evaluates the Adams-Bashforth predictor.
        procedure, public :: correct => adams_correct
            !! Applies the Adams-Moulton corrector.
        procedure, public :: evaluate_derivative => adams_evaluate_derivative
            !! Evaluates the derivative, resolving any mass matrix.
        procedure, public :: prepare_for_solve => adams_prepare_for_solve
            !! Discards the cached mass matrix factorization.
        procedure, public :: get_maximum_supported_order => adams_max_order
            !! Gets the highest order the method is permitted to use.
        procedure, public :: get_error_constant => adams_error_constant
            !! Gets the error constant for the requested order.
        procedure, public :: get_truncation_constant => adams_truncation_constant
            !! Gets the truncation error constant for the requested order.
        procedure, private :: initialize_workspace => adams_init_workspace
            !! Allocates the internal workspace arrays.
    end type

contains
! ------------------------------------------------------------------------------
pure function adams_max_order(this) result(rst)
    !! Gets the highest order this integrator is permitted to use.
    class(adams), intent(in) :: this
        !! The adams object.
    integer(int32) :: rst
        !! The order limit.
    rst = adams_order_limit
end function

! ------------------------------------------------------------------------------
pure function adams_error_constant(this, k) result(rst)
    !! Gets the constant relating the difference between the corrected and
    !! predicted solutions to the local truncation error of the formula.
    !!
    !! Writing the errors of the two formulae as \(C_p h^{k+1} y^{(k+1)}\) and
    !! \(C_c h^{k+1} y^{(k+1)}\), their difference eliminates the derivative
    !! and leaves the corrector error as
    !! \(\frac{C_c}{C_p - C_c}\left(y^c - y^p\right)\).  This is the Milne
    !! device, and the values returned here are those ratios.
    !!
    !! In terms of the Adams coefficients the ratio is
    !! \(\left|\gamma_k^*\right| / \gamma_{k-1}\), the denominator following
    !! from \(\gamma_k^* = \gamma_k - \gamma_{k-1}\).
    class(adams), intent(in) :: this
        !! The adams object.
    integer(int32), intent(in) :: k
        !! The order of interest.
    real(real64) :: rst
        !! The error constant.
    real(real64), parameter :: c(adams_order_limit) = [ &
        1.0d0 / 2.0d0, &
        1.0d0 / 6.0d0, &
        1.0d0 / 1.0d1, &
        1.9d1 / 2.7d2, &
        2.7d1 / 5.02d2, &
        8.63d2 / 1.995d4, &
        1.375d3 / 3.8174d4, &
        3.3953d4 / 1.10397d6, &
        5.7281d4 / 2.140034d6, &
        3.250433d6 / 1.37461698d8, &
        1.135053d6 / 5.3684506d7, &
        1.3695779093d10 / 7.1730003345d11 &
    ]
    rst = c(min(max(k, 1), adams_order_limit))
end function

! ------------------------------------------------------------------------------
pure function adams_truncation_constant(this, k) result(rst)
    !! Gets the constant relating the \((k+1)\)-th difference of the solution
    !! to the local truncation error of the corrector.
    !!
    !! These are the error constants of the Adams-Moulton formulae themselves,
    !! written \(\gamma_k^*\), and they differ from those returned by
    !! get_error_constant because the Adams predictor is not an extrapolation
    !! of the solution history.
    class(adams), intent(in) :: this
        !! The adams object.
    integer(int32), intent(in) :: k
        !! The order of interest.
    real(real64) :: rst
        !! The truncation error constant.
    real(real64), parameter :: c(adams_order_limit) = [ &
        1.0d0 / 2.0d0, &
        1.0d0 / 1.2d1, &
        1.0d0 / 2.4d1, &
        1.9d1 / 7.2d2, &
        3.0d0 / 1.6d2, &
        8.63d2 / 6.048d4, &
        2.75d2 / 2.4192d4, &
        3.3953d4 / 3.6288d6, &
        8.183d3 / 1.0368d6, &
        3.250433d6 / 4.790016d8, &
        4.671d3 / 7.8848d5, &
        1.3695779093d10 / 2.615348736d12 &
    ]
    rst = c(min(max(k, 1), adams_order_limit))
end function

! ------------------------------------------------------------------------------
subroutine adams_init_workspace(this, n)
    !! Allocates the internal workspace arrays.
    class(adams), intent(inout) :: this
        !! The adams object.
    integer(int32), intent(in) :: n
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%m_fpred)) then
        if (size(this%m_fpred) == n) return
        deallocate(this%m_fpred)
        deallocate(this%m_lu)
        deallocate(this%m_pivot)
        this%m_massFactored = .false.
    end if
    allocate(this%m_fpred(n), source = 0.0d0)
    allocate(this%m_lu(n, n), source = 0.0d0)
    allocate(this%m_pivot(n), source = 0)
end subroutine

! ------------------------------------------------------------------------------
subroutine adams_prepare_for_solve(this)
    !! Discards the cached mass matrix factorization.
    class(adams), intent(inout) :: this
        !! The adams object.
    this%m_massFactored = .false.
end subroutine

! ------------------------------------------------------------------------------
subroutine adams_evaluate_derivative(this, sys, x, y, dydx, args)
    !! Evaluates the derivative, resolving any mass matrix.
    !!
    !! An Adams method integrates the derivative itself, so a system written as
    !! \(M y' = f\) requires that the mass matrix be resolved rather than
    !! simply carried along.  The matrix is factored once and reused for as
    !! long as it remains valid, and a factorization that proves to be singular
    !! is reported rather than resolved in some least-squares sense, since a
    !! singular mass matrix places the problem beyond what this method can
    !! represent.
    class(adams), intent(inout) :: this
        !! The adams object.
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

    ! Local Variables
    integer(int32) :: i, n, flag
    real(real64) :: dmax, dmin

    ! Initialization
    n = size(y)
    call this%initialize_workspace(n)

    ! Evaluate the right-hand side
    call sys%fcn(x, y, dydx, args)
    if (.not.associated(sys%mass_matrix)) return

    ! Factor the mass matrix
    if (.not.this%m_massFactored .or. sys%get_is_mass_matrix_dependent()) then
        call sys%mass_matrix(x, y, this%m_lu, args)
        call DGETRF(n, n, this%m_lu, n, this%m_pivot, flag)
        if (flag /= 0) error stop DIFFEQ_SINGULAR_MATRIX_ERROR

        ! A vanishing pivot relative to the largest indicates that the matrix
        ! is singular to within the precision available
        dmax = 0.0d0
        dmin = huge(1.0d0)
        do i = 1, n
            dmax = max(dmax, abs(this%m_lu(i,i)))
            dmin = min(dmin, abs(this%m_lu(i,i)))
        end do
        if (dmin <= n * epsilon(1.0d0) * dmax) then
            error stop DIFFEQ_SINGULAR_MATRIX_ERROR
        end if
        this%m_massFactored = .true.
    end if

    ! Resolve the mass matrix
    call DGETRS('N', n, 1, this%m_lu, n, this%m_pivot, dydx, n, flag)
    if (flag /= 0) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
end subroutine

! ------------------------------------------------------------------------------
subroutine adams_predictor(this, sys, h, x, y, f, yn, fn, args)
    !! Evaluates the Adams-Bashforth predictor.
    !!
    !! The predictor integrates the polynomial fitted to the k most recent
    !! derivatives across the step.  The derivative is then evaluated at the
    !! predicted solution, which is the first evaluation of the PECE cycle and
    !! is retained for the corrector.
    class(adams), intent(inout) :: this
        !! The adams object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! The current value(s) of the dependent variable(s).
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives at x.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where this routine will write the next solution
        !! estimate at x + h.
    real(real64), intent(out), dimension(:) :: fn
        !! An N-element array where this routine will write the next
        !! derivative estimate at x + h.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information in and out
        !! of the differential equation subroutine.

    ! Local Variables
    integer(int32) :: i, k, n
    real(real64), allocatable, dimension(:) :: xv
    real(real64), allocatable, dimension(:,:) :: fv, dd

    ! Initialization
    n = size(y)
    k = max(1, min(this%get_order(), this%get_stored_step_count()))
    call this%initialize_workspace(n)

    ! The mass matrix is refreshed once per step attempt
    this%m_massFactored = .false.

    ! Gather the derivative history, most recent point first
    allocate(xv(k), fv(n, k), dd(n, k))
    do i = 1, k
        xv(i) = this%get_stored_variable(i)
        fv(:,i) = this%get_stored_derivative(i)
    end do

    ! Integrate the fitted polynomial across the step
    call newton_divided_differences(xv, fv, dd)
    call integrate_newton_polynomial(xv, dd, x, x + h, yn)
    yn = y + yn

    ! Evaluate the derivative at the prediction and retain it
    call this%evaluate_derivative(sys, x + h, yn, fn, args)
    this%m_fpred = fn
end subroutine

! ------------------------------------------------------------------------------
subroutine adams_correct(this, sys, h, x, ypred, ynew, converged, args)
    !! Applies the Adams-Moulton corrector.
    !!
    !! The corrector repeats the integration with the derivative at the new
    !! point included among the interpolation nodes, raising the degree of the
    !! fitted polynomial by one over the predictor at no additional cost, as
    !! that derivative was already evaluated.  It is applied exactly once; a
    !! step for which one correction does not suffice is rejected by the error
    !! test and retried on a shorter step.
    class(adams), intent(inout) :: this
        !! The adams object.
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
        !! True if the corrector succeeded; else, false.
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information in and out
        !! of the differential equation subroutine.

    ! Local Variables
    integer(int32) :: i, k, n
    real(real64), allocatable, dimension(:) :: xv
    real(real64), allocatable, dimension(:,:) :: fv, dd

    ! Initialization
    n = size(ypred)
    k = max(1, min(this%get_order(), this%get_stored_step_count()))

    ! The corrector interpolates the new point together with the k - 1 most
    ! recent stored points
    allocate(xv(k), fv(n, k), dd(n, k))
    xv(1) = x + h
    fv(:,1) = this%m_fpred
    do i = 2, k
        xv(i) = this%get_stored_variable(i - 1)
        fv(:,i) = this%get_stored_derivative(i - 1)
    end do

    ! Integrate the fitted polynomial across the step
    call newton_divided_differences(xv, fv, dd)
    call integrate_newton_polynomial(xv, dd, x, x + h, ynew)
    ynew = this%get_stored_solution(1) + ynew

    ! A correction that has ceased to be finite indicates that the step size
    ! has outrun the stability of the fixed-point iteration
    converged = all(abs(ynew) <= huge(1.0d0))
end subroutine

! ------------------------------------------------------------------------------
pure subroutine integrate_newton_polynomial(x, d, a, b, y)
    !! Integrates a polynomial stored in Newton divided-difference form.
    real(real64), intent(in), dimension(:) :: x
        !! An M-element array containing the interpolation nodes.
    real(real64), intent(in), dimension(:,:) :: d
        !! An N-by-M matrix containing the divided-difference coefficients.
    real(real64), intent(in) :: a
        !! The lower limit of integration.
    real(real64), intent(in) :: b
        !! The upper limit of integration.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where the result will be written.

    ! Gauss-Legendre quadrature is exact through degree fifteen, whereas the
    ! highest supported order produces a polynomial of degree eleven, so the
    ! integration introduces no error of its own.
    real(real64), parameter :: gp(8) = [ &
        -9.602898564975363d-1, -7.966664774136267d-1, &
        -5.255324099163290d-1, -1.834346424956498d-1, &
         1.834346424956498d-1,  5.255324099163290d-1, &
         7.966664774136267d-1,  9.602898564975363d-1 ]
    real(real64), parameter :: gw(8) = [ &
        1.012285362903763d-1, 2.223810344533745d-1, &
        3.137066458778873d-1, 3.626837833783620d-1, &
        3.626837833783620d-1, 3.137066458778873d-1, &
        2.223810344533745d-1, 1.012285362903763d-1 ]

    ! Local Variables
    integer(int32) :: i
    real(real64) :: c, hw
    real(real64), dimension(size(y)) :: p

    ! Process
    c = 0.5d0 * (a + b)
    hw = 0.5d0 * (b - a)
    y = 0.0d0
    do i = 1, size(gp)
        call evaluate_newton_polynomial(x, d, c + hw * gp(i), p)
        y = y + gw(i) * p
    end do
    y = hw * y
end subroutine

! ------------------------------------------------------------------------------
end module
