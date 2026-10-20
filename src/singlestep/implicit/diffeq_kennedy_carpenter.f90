! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module diffeq_kennedy_carpenter
    use iso_fortran_env
    use diffeq_base
    use diffeq_singlestep
    use diffeq_errors
    use linalg
    implicit none
    private
    public :: kennedy_carpenter
    public :: kennedy_carpenter_4
    public :: kennedy_carpenter_5

    real(real64), parameter :: kc_jacobian_refresh_rate = 3.0d-2
        ! A Newton contraction rate above this value marks the Jacobian as
        ! stale, so it is re-evaluated before the next step.

    type, abstract, extends(single_step_integrator) :: kennedy_carpenter
        !! Shared implementation for Kennedy--Carpenter ESDIRK methods.
        !!
        !! The stage equations are solved with a simplified Newton iteration.
        !! A supplied mass matrix is incorporated in each stage system as
        !! \(M-h \gamma J\), while the embedded tableau supplies the error
        !! estimate used by the inherited adaptive driver.  These solvers use
        !! the inherited configurable PI step-size controller.
        !!
        !! Every implicit stage shares the diagonal coefficient \(\gamma\), so
        !! a single factorization of \(M-h \gamma J\) serves every stage and
        !! Newton iteration within a step.  The factorization is retained
        !! until the step size changes, and the Jacobian is retained across
        !! steps until the Newton iteration converges slowly or fails.  A
        !! failed iteration first refreshes any stale Jacobian or
        !! factorization; if the iteration still fails, the step is rejected
        !! and retried with a smaller step size.
        !!
        !! The Newton iteration is governed by a user-defined iteration limit
        !! and convergence tolerance.  The applied tolerance never exceeds one
        !! percent of the smallest positive integration tolerance.
        real(real64), private, allocatable, dimension(:,:) :: m_mass
        real(real64), private, allocatable, dimension(:) :: m_startDerivative
        logical, private :: m_massComputed = .false.
        integer, private :: m_newtonIterLimit = 12
            ! The limit on the number of Newton iterations per stage.
        real(real64), private :: m_newtonTolerance = 1.0d-8
            ! The convergence tolerance for the Newton iteration.
        real(real64), private, allocatable, dimension(:,:) :: m_jac
            ! The Jacobian, less any mass-matrix state derivative, used in
            ! the iteration matrix.
        real(real64), private, allocatable, dimension(:,:) :: m_qr
            ! The QR factored iteration matrix.
        real(real64), private, allocatable, dimension(:) :: m_tau
            ! The QR factorization scalar factors.
        integer, private, allocatable, dimension(:) :: m_pivot
            ! The QR factorization column pivots.
        real(real64), private :: m_factorStep = 0.0d0
            ! The value of h * gamma used to form the factored matrix.
        logical, private :: m_haveFactor = .false.
            ! True if a usable factorization is stored.
        logical, private :: m_jacCurrent = .false.
            ! True if the Jacobian was evaluated at the current step start.
        logical, private :: m_refreshJac = .true.
            ! True if the Jacobian must be re-evaluated before the next step.
    contains
        procedure, public :: solve => kc_solve
        procedure, public :: get_newton_step_limit => kc_get_newton_step_limit
            !! Gets the maximum number of Newton iterations allowed per stage.
        procedure, public :: set_newton_step_limit => kc_set_newton_step_limit
            !! Sets the maximum number of Newton iterations allowed per stage.
        procedure, public :: get_newton_tolerance => kc_get_newton_tolerance
            !! Gets the Newton iteration convergence tolerance.
        procedure, public :: set_newton_tolerance => kc_set_newton_tolerance
            !! Sets the Newton iteration convergence tolerance.
        procedure, public :: pre_step_action => kc_pre_step
            !! Performs any pre-step actions.
        procedure, public :: attempt_step => kc_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => kc_post_step
            !! Performs any post-step actions.
        procedure, public :: interpolate => kc_interpolate
            !! Performs the interpolation.
        procedure, public :: get_is_fsal => kc_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => kc_get_stage_count
            !! Gets the stage count for this integrator.
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_4
        !! Fourth-order ARK4(3)6L[2]SA ESDIRK method.
        !! The inherited PI step-size controller uses the embedded error
        !! estimate to select the next step size.
    contains
        procedure, public :: get_order => kc4_get_order
            !! Gets the order of the integrator.
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_5
        !! Fifth-order ARK5(4)8L[2]SA ESDIRK method.
        !! The inherited PI step-size controller uses the embedded error
        !! estimate to select the next step size.
    contains
        procedure, public :: get_order => kc5_get_order
            !! Gets the order of the integrator.
    end type

contains
! ------------------------------------------------------------------------------
subroutine kc_solve(this, sys, x, iv, args)
    class(kennedy_carpenter), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(*), intent(inout), optional, target :: args

    this%m_massComputed = .false.
    this%m_haveFactor = .false.
    this%m_jacCurrent = .false.
    this%m_refreshJac = .true.
    this%m_factorStep = 0.0d0
    if (allocated(this%m_mass)) deallocate(this%m_mass)
    if (allocated(this%m_startDerivative)) deallocate(this%m_startDerivative)
    call ssi_ode_solver(this, sys, x, iv, args)
end subroutine

! ------------------------------------------------------------------------------
pure function kc_get_newton_step_limit(this) result(rst)
    !! Gets the maximum number of Newton iterations allowed per stage.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    integer :: rst
        !! The iteration limit.
    rst = this%m_newtonIterLimit
end function

! --------------------
subroutine kc_set_newton_step_limit(this, x)
    !! Sets the maximum number of Newton iterations allowed per stage.  A
    !! value less than one causes DIFFEQ_INVALID_INPUT_ERROR.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    integer, intent(in) :: x
        !! The iteration limit.
    if (x < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_newtonIterLimit = x
end subroutine

! ------------------------------------------------------------------------------
pure function kc_get_newton_tolerance(this) result(rst)
    !! Gets the convergence tolerance on the estimated Newton iteration error,
    !! relative to the larger of one and the norm of the stage state.  The
    !! tolerance applied during integration is the lesser of this value and
    !! one percent of the smallest positive integration tolerance.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    real(real64) :: rst
        !! The convergence tolerance.
    rst = this%m_newtonTolerance
end function

! --------------------
subroutine kc_set_newton_tolerance(this, x)
    !! Sets the convergence tolerance on the estimated Newton iteration error,
    !! relative to the larger of one and the norm of the stage state.  The
    !! tolerance applied during integration is the lesser of this value and
    !! one percent of the smallest positive integration tolerance.  A value
    !! that is not positive and finite causes DIFFEQ_INVALID_INPUT_ERROR.
    use ieee_arithmetic, only : ieee_is_finite
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(in) :: x
        !! The convergence tolerance.
    if (.not.ieee_is_finite(x)) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (x <= 0.0d0) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_newtonTolerance = x
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_pre_step(this, prevs, sys, h, x, y, f, args)
    !! Marks the stored Jacobian as belonging to an earlier step once the
    !! solution has advanced.
    !!
    !! The Jacobian itself is formed lazily by the step attempt, as it is only
    !! re-evaluated when the Newton iteration requires it.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
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
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Process
    if (prevs) this%m_jacCurrent = .false.
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_post_step(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    !! Retains the first-stage derivative for cubic Hermite interpolation.
    !!
    !! For a mass-matrix problem, the driver's f is a raw RHS rather than a
    !! derivative. The first stage supplies the consistent left-end slope;
    !! the stiffly accurate final stage already supplies the right-end slope.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
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
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    this%m_startDerivative = k(:,1)
end subroutine

! ------------------------------------------------------------------------------
pure function kc_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .false.
end function

! ------------------------------------------------------------------------------
pure function kc_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    integer :: rst
        !! The stage count.
    select type (this)
    type is (kennedy_carpenter_4)
        rst = 6
    type is (kennedy_carpenter_5)
        rst = 8
    class default
        rst = 0
    end select
end function

! ------------------------------------------------------------------------------
pure function kc4_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(kennedy_carpenter_4), intent(in) :: this
        !! The kennedy_carpenter_4 object.
    integer :: rst
        !! The order.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function kc5_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(kennedy_carpenter_5), intent(in) :: this
        !! The kennedy_carpenter_5 object.
    integer :: rst
        !! The order.
    rst = 5
end function

! ------------------------------------------------------------------------------
subroutine kc_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    !! Attempts an integration step for this integrator.
    !!
    !! The state at each stage is solved from
    !! \[
    !! M(Y_i-Y_i^*)-h a_{ii}f(x+c_i h,Y_i)=0,
    !! \]
    !! using Newton iteration.  The high- and embedded-order solutions are
    !! \(y_{n+1}=y_n+h\sum b_i k_i\) and
    !! \(\widehat y_{n+1}=y_n+h\sum d_i k_i\), respectively.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
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
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: i, j, n, stages
    logical :: usemass, massDependent, converged
    real(real64) :: a(8,8), b(8), d(8), c(8), hg, theta, thetaMax
    real(real64), allocatable, dimension(:) :: base, embedded

    ! Initialization
    n = size(y)
    usemass = sys%has_mass_matrix()
    massDependent = sys%get_is_mass_matrix_dependent()
    call kc_table(this, a, b, d, c, stages)
    hg = h * a(2,2)
    call kc_init_workspace(this, n)
    allocate(base(n), embedded(n))

    ! Process
    ! The first stage of an ESDIRK method is explicit
    k = 0.0d0
    if (usemass) then
        if (allocated(this%m_mass)) then
            if (size(this%m_mass, 1) /= n) then
                deallocate(this%m_mass)
                this%m_massComputed = .false.
            end if
        end if
        if (.not.allocated(this%m_mass)) allocate(this%m_mass(n, n))
        if (.not.this%m_massComputed .or. massDependent) then
            call sys%compute_mass_matrix(x, y, this%m_mass, args)
            this%m_massComputed = .true.
            if (massDependent) this%m_haveFactor = .false.
        end if
        call kc_consistent_derivative(this, sys, x, y, f, massDependent, &
            k(:,1), args)
    else
        k(:,1) = f
    end if

    ! Every implicit stage shares the diagonal coefficient, so one factored
    ! iteration matrix serves the entire step.  A Newton failure computed with
    ! a Jacobian from an earlier step refreshes the Jacobian and retries; a
    ! failure with a current Jacobian rejects the step.
    do
        if (this%m_refreshJac .and. .not.this%m_jacCurrent) then
            call kc_update_jacobian(this, sys, x, y, k(:,1), usemass, &
                massDependent, args)
        end if
        if (.not.this%m_haveFactor .or. this%m_factorStep /= hg) then
            call kc_factor(this, hg, usemass)
        end if

        thetaMax = 0.0d0
        converged = .true.
        do i = 2, stages
            ! Accumulate the contributions from the previously computed stages
            base = y
            do j = 1, i - 1
                base = base + h * a(i,j) * k(:,j)
            end do

            call kc_newton_stage(this, sys, x + c(i) * h, base, hg, &
                k(:,i-1), usemass, massDependent, k(:,i), converged, &
                theta, args)
            thetaMax = max(thetaMax, theta)
            if (.not.converged) exit
        end do
        if (converged) exit

        if (this%m_jacCurrent) then
            call kc_reject_step(this, y, yn, yerr)
            fn = k(:,1)
            return
        end if
        this%m_refreshJac = .true.
    end do

    ! Slow convergence indicates the Jacobian no longer represents the
    ! problem well, so refresh it before the next step
    if (thetaMax > kc_jacobian_refresh_rate) this%m_refreshJac = .true.

    ! Form the solution estimate and the embedded solution estimate
    yn = y
    embedded = y
    do i = 1, stages
        yn = yn + h * b(i) * k(:,i)
        embedded = embedded + h * d(i) * k(:,i)
    end do
    yerr = yn - embedded

    ! Both tableaus are stiffly accurate.  The final stage derivative is
    ! therefore the derivative at the accepted state, even when M is singular.
    fn = k(:, stages)
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_init_workspace(this, n)
    !! Allocates the Jacobian and factorization workspace.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    integer, intent(in) :: n
        !! The number of equations.

    if (allocated(this%m_jac)) then
        if (size(this%m_jac, 1) == n) return
        deallocate(this%m_jac, this%m_qr, this%m_tau, this%m_pivot)
    end if
    allocate(this%m_jac(n,n), this%m_qr(n,n), this%m_tau(n), this%m_pivot(n))
    this%m_haveFactor = .false.
    this%m_jacCurrent = .false.
    this%m_refreshJac = .true.
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_update_jacobian(this, sys, x, y, derivative, usemass, &
    massDependent, args)
    !! Evaluates the Jacobian at the start of the current step.
    !!
    !! For a state-dependent mass matrix, the linearization of
    !! \(M(Y)(Y-Y^*)\) contributes \(\partial (M y')/\partial y\), which is
    !! approximated with the step-start derivative and subtracted from the
    !! stored Jacobian.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x
        !! The value of the independent variable at the start of the step.
    real(real64), intent(in), dimension(:) :: y
        !! The solution at the start of the step.
    real(real64), intent(in), dimension(:) :: derivative
        !! The derivative at the start of the step.
    logical, intent(in) :: usemass, massDependent
        !! Mass-matrix configuration flags.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: i, n
    real(real64) :: fdstep
    real(real64), allocatable :: mass_perturbed(:,:), state(:)

    n = size(y)
    call sys%compute_jacobian(x, y, this%m_jac, args)
    call this%log_jacobian_evaluation()
    if (sys%use_finite_difference_jacobian()) then
        call this%log_function_evaluation(n + 1)
    end if

    if (usemass .and. massDependent) then
        fdstep = sys%get_finite_difference_step()
        allocate(mass_perturbed(n,n), state(n))
        do i = 1, n
            state = y
            state(i) = state(i) + fdstep
            call sys%compute_mass_matrix(x, state, mass_perturbed, args)
            this%m_jac(:,i) = this%m_jac(:,i) - &
                matmul((mass_perturbed - this%m_mass) / fdstep, derivative)
        end do
    end if

    this%m_jacCurrent = .true.
    this%m_refreshJac = .false.
    this%m_haveFactor = .false.
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_factor(this, hg, usemass)
    !! Forms and factors the iteration matrix \(M - h \gamma J\).
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(in) :: hg
        !! The product of the step size and the diagonal coefficient.
    logical, intent(in) :: usemass
        !! True if a mass matrix is used; else, false.

    ! Local Variables
    integer :: i, n
    real(real64), allocatable :: system(:,:)

    n = size(this%m_jac, 1)
    if (usemass) then
        system = this%m_mass - hg * this%m_jac
    else
        system = -hg * this%m_jac
        do i = 1, n
            system(i,i) = system(i,i) + 1.0d0
        end do
    end if
    this%m_pivot = 0
    call qr_factor(system, this%m_pivot, tau = this%m_tau, qr = this%m_qr)
    call this%log_matrix_factorization()
    this%m_factorStep = hg
    this%m_haveFactor = .true.
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_newton_stage(this, sys, x, base, hg, guess, usemass, &
    massDependent, derivative, converged, thetaMax, args)
    !! Solves one stage equation
    !! \(M(Y_i-Y_i^*)-h \gamma f(x+c_i h,Y_i)=0\)
    !! by simplified Newton iteration using the stored factorization.
    !!
    !! The iteration converges linearly at a rate \(\theta\), so the error
    !! remaining after an update \(\Delta\) is estimated as
    !! \(\eta\|\Delta\|\), with \(\eta=\theta/(1-\theta)\).  The iteration
    !! has converged once this estimate falls below the Newton tolerance,
    !! relative to the larger of one and the norm of the stage state.  No
    !! contraction rate is available after the first iteration, so the update
    !! itself must then satisfy the tolerance.  The Newton tolerance is capped
    !! at one percent of the smallest positive integration tolerance.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x
        !! The value of the independent variable at the stage.
    real(real64), intent(in), dimension(:) :: base
        !! The explicit portion of the stage state, \(Y_i^*\).
    real(real64), intent(in) :: hg
        !! The product of the step size and the diagonal coefficient.
    real(real64), intent(in), dimension(:) :: guess
        !! The derivative used to predict the stage state.
    logical, intent(in) :: usemass, massDependent
        !! Mass-matrix configuration flags.
    real(real64), intent(out), dimension(:) :: derivative
        !! The stage derivative, \((Y_i-Y_i^*)/(h \gamma)\).
    logical, intent(out) :: converged
        !! True if the iteration converged; else, false.
    real(real64), intent(out) :: thetaMax
        !! The largest observed contraction rate of the iteration.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: iteration
    real(real64) :: tol, dnorm, dold, theta, eta
    real(real64), dimension(size(base)) :: state, f, residual, delta
    real(real64), allocatable :: mass(:,:)

    tol = kc_effective_newton_tolerance(this, size(base))
    if (massDependent) allocate(mass(size(base), size(base)))
    state = base + hg * guess
    converged = .false.
    thetaMax = 0.0d0
    dold = 0.0d0
    eta = 1.0d0
    do iteration = 1, this%get_newton_step_limit()
        call sys%compute_rhs(x, state, f, args)
        call this%log_function_evaluation()
        if (.not.usemass) then
            residual = state - base - hg * f
        else if (massDependent) then
            call sys%compute_mass_matrix(x, state, mass, args)
            residual = matmul(mass, state - base) - hg * f
        else
            residual = matmul(this%m_mass, state - base) - hg * f
        end if
        delta = solve_qr(this%m_qr, this%m_tau, this%m_pivot, -residual)
        state = state + delta
        dnorm = norm2(delta)
        if (iteration > 1) then
            theta = dnorm / dold
            thetaMax = max(thetaMax, theta)
            if (theta >= 1.0d0) exit
            eta = theta / (1.0d0 - theta)
        end if
        if (eta * dnorm <= tol * max(1.0d0, norm2(state))) then
            converged = .true.
            exit
        end if
        dold = dnorm
    end do
    derivative = (state - base) / hg
end subroutine

! ------------------------------------------------------------------------------
function kc_effective_newton_tolerance(this, n) result(rst)
    !! Gets the Newton tolerance actually applied, which is the lesser of the
    !! user-defined Newton tolerance and one percent of the smallest positive
    !! integration tolerance.  Iteration error larger than the integration
    !! tolerance would otherwise leave algebraic constraints inconsistent by
    !! an amount no step size reduction can remove.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    integer, intent(in) :: n
        !! The number of equations.
    real(real64) :: rst
        !! The tolerance.

    ! Local Variables
    integer :: i
    real(real64) :: atol, rtol

    rst = this%get_newton_tolerance()
    do i = 1, n
        atol = this%get_absolute_tolerance(i)
        rtol = this%get_relative_tolerance(i)
        if (atol > 0.0d0) rst = min(rst, 1.0d-2 * atol)
        if (rtol > 0.0d0) rst = min(rst, 1.0d-2 * rtol)
    end do
end function

! ------------------------------------------------------------------------------
subroutine kc_reject_step(this, y, yn, yerr)
    !! Produces an error estimate that forces the driver to reject the step
    !! and reduce the step size by roughly one half, following a Newton
    !! iteration failure.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(in), dimension(:) :: y
        !! The solution at the start of the step.
    real(real64), intent(out), dimension(:) :: yn
        !! The solution estimate, set to y.
    real(real64), intent(out), dimension(:) :: yerr
        !! The error estimate.

    ! Local Variables
    integer :: i
    real(real64) :: alpha, e

    ! The inherited controller scales the step by fs * e**(-alpha)
    alpha = 1.0d0 / this%get_order() - &
        0.75d0 * this%get_step_size_control_parameter()
    if (alpha > 0.0d0) then
        e = min(max(2.0d0, &
            (2.0d0 * this%get_step_size_factor())**(1.0d0 / alpha)), 1.0d10)
    else
        e = 1.0d4
    end if
    yn = y
    do i = 1, size(y)
        yerr(i) = e * (this%get_absolute_tolerance(i) + &
            this%get_relative_tolerance(i) * abs(y(i)))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_consistent_derivative(this, sys, x, y, f, massDependent, &
    derivative, args)
    !! Computes a consistent initial derivative for an index-1 DAE.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x, y(:), f(:)
        !! The step-start point and right-hand side.
    logical, intent(in) :: massDependent
        !! True if the mass matrix depends upon the state.
    real(real64), intent(out) :: derivative(:)
        !! The consistent derivative.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    real(real64), allocatable :: singular_values(:), left_vectors(:,:), &
        jac(:,:), fplus(:), augmented(:,:), rhs(:)
    real(real64) :: fdstep, rank_tolerance, scale
    integer :: n, rank, nullity, i

    n = size(y)
    call svd(this%m_mass, singular_values, left_vectors)
    scale = max(1.0d0, singular_values(1))
    rank_tolerance = 100.0d0 * epsilon(1.0d0) * scale
    rank = count(singular_values > rank_tolerance)
    call this%log_matrix_factorization()

    if (rank == n) then
        derivative = solve_kc_system(this%m_mass, f)
        return
    end if

    ! The algebraic rows require the true Jacobian at the step start.  For a
    ! constant mass matrix, the stored Jacobian is that matrix, so it is
    ! evaluated here once and shared with the stage iteration.
    nullity = n - rank
    allocate(jac(n,n), fplus(n), augmented(n,n), rhs(n))
    if (massDependent) then
        call sys%compute_jacobian(x, y, jac, args)
        call this%log_jacobian_evaluation()
        if (sys%use_finite_difference_jacobian()) then
            call this%log_function_evaluation(n + 1)
        end if
    else
        if (.not.this%m_jacCurrent) then
            call kc_update_jacobian(this, sys, x, y, f, .true., .false., args)
        end if
        jac = this%m_jac
    end if
    fdstep = sys%get_finite_difference_step()
    call sys%compute_rhs(x + fdstep, y, fplus, args)
    call this%log_function_evaluation()
    fplus = (fplus - f) / fdstep

    augmented = 0.0d0
    rhs = 0.0d0
    augmented(1:rank,:) = matmul(transpose(left_vectors(:,1:rank)), this%m_mass)
    rhs(1:rank) = matmul(transpose(left_vectors(:,1:rank)), f)
    do i = 1, nullity
        augmented(rank+i,:) = matmul(left_vectors(:,rank+i), jac)
        rhs(rank+i) = -dot_product(left_vectors(:,rank+i), fplus)
    end do

    derivative = solve_kc_system(augmented, rhs)
end subroutine

subroutine kc_table(this, a, b, d, c, stages)
    !! Populates the Butcher tableau for the requested integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(out) :: a(8,8)
        !! The matrix of stage coefficients.  Only the leading
        !! stages-by-stages block is populated.
    real(real64), intent(out) :: b(8)
        !! The array of weights defining the higher-order solution.
    real(real64), intent(out) :: d(8)
        !! The array of weights defining the embedded solution.
    real(real64), intent(out) :: c(8)
        !! The array of nodes at which each stage is evaluated.
    integer, intent(out) :: stages
        !! The number of stages used by the integrator.

    ! Initialization
    a = 0.0d0
    b = 0.0d0
    d = 0.0d0
    c = 0.0d0

    ! Process
    select type (this)
    type is (kennedy_carpenter_4)
        ! ARK4(3)6L[2]SA
        stages = 6

        a(2,1) = 1.0d0 / 4.0d0
        a(2,2) = 1.0d0 / 4.0d0

        a(3,1) = 8611.0d0 / 62500.0d0
        a(3,2) = -1743.0d0 / 31250.0d0
        a(3,3) = 1.0d0 / 4.0d0

        a(4,1) = 5012029.0d0 / 34652500.0d0
        a(4,2) = -654441.0d0 / 2922500.0d0
        a(4,3) = 174375.0d0 / 388108.0d0
        a(4,4) = 1.0d0 / 4.0d0

        a(5,1) = 15267082809.0d0 / 155376265600.0d0
        a(5,2) = -71443401.0d0 / 120774400.0d0
        a(5,3) = 730878875.0d0 / 902184768.0d0
        a(5,4) = 2285395.0d0 / 8070912.0d0
        a(5,5) = 1.0d0 / 4.0d0

        a(6,1) = 82889.0d0 / 524892.0d0
        a(6,3) = 15625.0d0 / 83664.0d0
        a(6,4) = 69875.0d0 / 102672.0d0
        a(6,5) = -2260.0d0 / 8211.0d0
        a(6,6) = 1.0d0 / 4.0d0

        ! The method is stiffly accurate, so the weights match the last stage
        b(1) = a(6,1)
        b(3) = a(6,3)
        b(4) = a(6,4)
        b(5) = a(6,5)
        b(6) = a(6,6)

        d(1) = 4586570599.0d0 / 29645900160.0d0
        d(3) = 178811875.0d0 / 945068544.0d0
        d(4) = 814220225.0d0 / 1159782912.0d0
        d(5) = -3700637.0d0 / 11593932.0d0
        d(6) = 61727.0d0 / 225920.0d0

        c(2) = 0.5d0
        c(3) = 83.0d0 / 250.0d0
        c(4) = 31.0d0 / 50.0d0
        c(5) = 17.0d0 / 20.0d0
        c(6) = 1.0d0
    type is (kennedy_carpenter_5)
        ! ARK5(4)8L[2]SA
        stages = 8

        a(2,1) = 41.0d0 / 200.0d0
        a(2,2) = 41.0d0 / 200.0d0

        a(3,1) = 41.0d0 / 400.0d0
        a(3,2) = -567603406766.0d0 / 11931857230679.0d0
        a(3,3) = 41.0d0 / 200.0d0

        a(4,1) = 683785636431.0d0 / 9252920307686.0d0
        a(4,3) = -110385047103.0d0 / 1367015193373.0d0
        a(4,4) = 41.0d0 / 200.0d0

        a(5,1) = 3016520224154.0d0 / 10081342136671.0d0
        a(5,3) = 30586259806659.0d0 / 12414158314087.0d0
        a(5,4) = -22760509404356.0d0 / 11113319521817.0d0
        a(5,5) = 41.0d0 / 200.0d0

        a(6,1) = 218866479029.0d0 / 1489978393911.0d0
        a(6,3) = 638256894668.0d0 / 5436446318841.0d0
        a(6,4) = -1179710474555.0d0 / 5321154724896.0d0
        a(6,5) = -60928119172.0d0 / 8023461067671.0d0
        a(6,6) = 41.0d0 / 200.0d0

        a(7,1) = 1020004230633.0d0 / 5715676835656.0d0
        a(7,3) = 25762820946817.0d0 / 25263940353407.0d0
        a(7,4) = -2161375909145.0d0 / 9755907335909.0d0
        a(7,5) = -211217309593.0d0 / 5846859502534.0d0
        a(7,6) = -4269925059573.0d0 / 7827059040749.0d0
        a(7,7) = 41.0d0 / 200.0d0

        a(8,1) = -872700587467.0d0 / 9133579230613.0d0
        a(8,4) = 22348218063261.0d0 / 9555858737531.0d0
        a(8,5) = -1143369518992.0d0 / 8141816002931.0d0
        a(8,6) = -39379526789629.0d0 / 19018526304540.0d0
        a(8,7) = 32727382324388.0d0 / 42900044865799.0d0
        a(8,8) = 41.0d0 / 200.0d0

        ! The method is stiffly accurate, so the weights match the last stage
        b(1) = a(8,1)
        b(4) = a(8,4)
        b(5) = a(8,5)
        b(6) = a(8,6)
        b(7) = a(8,7)
        b(8) = a(8,8)

        d(1) = -975461918565.0d0 / 9796059967033.0d0
        d(4) = 78070527104295.0d0 / 32432590147079.0d0
        d(5) = -548382580838.0d0 / 3424219808633.0d0
        d(6) = -33438840321285.0d0 / 15594753105479.0d0
        d(7) = 3629800801594.0d0 / 4656183773603.0d0
        d(8) = 4035322873751.0d0 / 18575991585200.0d0

        c(2) = 41.0d0 / 100.0d0
        c(3) = 2935347310677.0d0 / 11292855782101.0d0
        c(4) = 1426016391358.0d0 / 7196633302097.0d0
        c(5) = 0.92d0
        c(6) = 0.24d0
        c(7) = 0.6d0
        c(8) = 1.0d0
    end select
end subroutine

! ------------------------------------------------------------------------------
function solve_kc_system(a, b) result(x)
    !! Solves the linear system \( A x = b \) by means of a QR factorization
    !! with column pivoting.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N system matrix.
    real(real64), intent(in), dimension(:) :: b
        !! An N-element array containing the right-hand side.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the solution.

    ! Process
    x = solve_least_squares_full(a, b)
end function

! ------------------------------------------------------------------------------
subroutine kc_interpolate(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    !!
    !! A cubic Hermite polynomial is constructed from the solution and
    !! derivative values at each end of the step.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
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
    real(real64) :: h, s

    ! Initialization
    h = xn1 - xn
    s = (x - xn) / h

    ! Process
    y = (2.0d0 * s**3 - 3.0d0 * s**2 + 1.0d0) * yn + &
        (-2.0d0 * s**3 + 3.0d0 * s**2) * yn1 + &
        (s**3 - s**2) * h * fn1
    if (allocated(this%m_startDerivative)) then
        y = y + (s**3 - 2.0d0 * s**2 + s) * h * this%m_startDerivative
    else
        y = y + (s**3 - 2.0d0 * s**2 + s) * h * fn
    end if
end subroutine

! ------------------------------------------------------------------------------
end module