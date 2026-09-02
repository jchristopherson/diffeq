module diffeq_bdf
    !! Implicit BDF solvers.
    !!
    !! Backward differentiation formulae replace the derivative in
    !! \[
    !! M(x,y) y' = f(x,y)
    !! \]
    !! with the derivative, evaluated at \(x_{n+1}\), of the polynomial that
    !! interpolates the new solution point together with the \(k\) points
    !! preceding it.  The resulting formula
    !! \[
    !! M(x_{n+1},y_{n+1}) \sum_{j=0}^{k} c_j y_{n+1-j} = f(x_{n+1},y_{n+1})
    !! \]
    !! is implicit in \(y_{n+1}\) and is solved by Newton iteration.  Because
    !! the coefficients are derived from the actual spacing of the stored
    !! points, a change of step size is handled exactly rather than by
    !! interpolating the history onto a uniform grid.
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use diffeq_multistep
    implicit none
    private
    public :: bdf

    integer(int32), parameter :: bdf_order_limit = 5
        ! BDF formulae are only zero-stable through order six, and the sixth
        ! order formula has so small a region of stability as to be of little
        ! practical use.  The order is therefore capped at five.

    type, extends(multistep_integrator) :: bdf
        !! A variable step-size, variable order (1-5) BDF integrator for stiff
        !! systems and index-1 differential-algebraic equations.
        !!
        !! The iteration matrix formed for the Newton corrector is
        !! \[
        !! A = \frac{\alpha_0}{h_n} M - J, \qquad
        !! J = \frac{\partial f}{\partial y},
        !! \]
        !! which is assembled and factored rather than requiring \(M\) to be
        !! inverted.  A singular mass matrix is therefore permitted: rows of
        !! \(M\) that are identically zero impose algebraic constraints, and
        !! the problem becomes a differential-algebraic equation.  So long as
        !! that system is of index one, \(A\) remains nonsingular and the
        !! solver proceeds normally.
        !!
        !! The method restarts at first order, so it does not require a
        !! consistent initial derivative.  Initial conditions that satisfy the
        !! algebraic constraints are still recommended, as the first step must
        !! otherwise absorb the inconsistency.
        real(real64), private :: m_c0 = 0.0d0
            ! The leading BDF coefficient, alpha_0 / h, for the current step.
        real(real64), private, allocatable, dimension(:) :: m_histTerm
            ! The contribution of the stored solution history to the formula.
        real(real64), private, allocatable, dimension(:,:) :: m_jac
            ! The Jacobian matrix.
        real(real64), private, allocatable, dimension(:,:) :: m_mass
            ! The mass matrix.
        real(real64), private, allocatable, dimension(:) :: m_rhs
            ! The value of f(x,y).
        real(real64), private, allocatable, dimension(:) :: m_ydot
            ! The BDF approximation of the derivative at x + h.
        logical, private :: m_massComputed = .false.
            ! True if the mass matrix has been computed for this step.
    contains
        procedure, public :: predictor => bdf_predictor
            !! Evaluates the predictor and forms the coefficients of the step.
        procedure, public :: build_linear_system => bdf_build_linear_system
            !! Builds the linearized system for the implicit corrector.
        procedure, public :: get_maximum_supported_order => bdf_max_order
            !! Gets the highest order the method is permitted to use.
        procedure, public :: get_error_constant => bdf_error_constant
            !! Gets the error constant for the requested order.
        procedure, private :: initialize_workspace => bdf_init_workspace
            !! Allocates the internal workspace arrays.
    end type

contains
! ------------------------------------------------------------------------------
pure function bdf_max_order(this) result(rst)
    !! Gets the highest order this integrator is permitted to use.
    class(bdf), intent(in) :: this
        !! The bdf object.
    integer(int32) :: rst
        !! The order limit.
    rst = bdf_order_limit
end function

! ------------------------------------------------------------------------------
pure function bdf_error_constant(this, k) result(rst)
    !! Gets the constant relating the difference between the corrected and
    !! predicted solutions to the local truncation error of the formula.
    !!
    !! The predictor of order \(k\) and the corrector of order \(k\) differ by
    !! a quantity proportional to the leading error term of the formula, and
    !! for a BDF method that constant of proportionality is \(1/(k+1)\).
    class(bdf), intent(in) :: this
        !! The bdf object.
    integer(int32), intent(in) :: k
        !! The order of interest.
    real(real64) :: rst
        !! The error constant.
    rst = 1.0d0 / (k + 1.0d0)
end function

! ------------------------------------------------------------------------------
subroutine bdf_init_workspace(this, n)
    !! Allocates the internal workspace arrays.
    class(bdf), intent(inout) :: this
        !! The bdf object.
    integer(int32), intent(in) :: n
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%m_jac)) then
        if (size(this%m_jac, 1) == n) return
        deallocate(this%m_jac)
        deallocate(this%m_mass)
        deallocate(this%m_rhs)
        deallocate(this%m_ydot)
        deallocate(this%m_histTerm)
    end if
    allocate( &
        this%m_jac(n, n), &
        this%m_mass(n, n), &
        this%m_rhs(n), &
        this%m_ydot(n), &
        this%m_histTerm(n) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine bdf_predictor(this, sys, h, x, y, f, yn, fn, args)
    !! Evaluates the predictor and forms the coefficients for the step.
    !!
    !! The prediction is the extrapolation to \(x + h\) of the polynomial
    !! passing through the stored solution history.  This routine also builds
    !! the BDF coefficients for the step and collapses the history into the
    !! single vector the corrector requires, both of which are fixed for the
    !! duration of the step and so are computed only once here rather than on
    !! every Newton iteration.
    class(bdf), intent(inout) :: this
        !! The bdf object.
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
    integer(int32) :: i, k, n, nc, np
    real(real64), allocatable, dimension(:) :: xv, c, nodes
    real(real64), allocatable, dimension(:,:) :: yv, dd

    ! Initialization
    n = size(y)
    nc = this%get_stored_step_count()
    k = max(1, min(this%get_order(), nc))
    call this%initialize_workspace(n)

    ! Gather the solution history, most recent point first
    np = min(k + 1, nc)
    allocate(xv(np), yv(n, np), dd(n, np))
    do i = 1, np
        xv(i) = this%get_stored_variable(i)
        yv(:,i) = this%get_stored_solution(i)
    end do

    ! Extrapolate the history to the end of the step.  On the very first step
    ! only the initial condition is available, so the prediction reduces to
    ! the current solution.
    if (np == 1) then
        yn = yv(:,1)
    else
        call newton_divided_differences(xv, yv, dd)
        call evaluate_newton_polynomial(xv, dd, x + h, yn)
    end if
    call sys%fcn(x + h, yn, fn, args)

    ! Form the coefficients of the formula from the actual node locations
    allocate(nodes(k + 1), c(k + 1))
    nodes(1) = x + h
    nodes(2:k+1) = xv(1:k)
    call bdf_coefficients(nodes, c)

    ! Collapse the history into the term the corrector residual requires
    this%m_c0 = c(1)
    this%m_histTerm = 0.0d0
    do i = 1, k
        this%m_histTerm = this%m_histTerm + c(i+1) * yv(:,i)
    end do
    this%m_massComputed = .false.
end subroutine

! ------------------------------------------------------------------------------
pure subroutine bdf_coefficients(t, c)
    !! Computes the coefficients of the BDF formula for the supplied nodes.
    !!
    !! The coefficients are the derivatives, evaluated at the newest node, of
    !! the Lagrange basis polynomials belonging to the supplied nodes, so that
    !! \( P'(t_1) = \sum_j c_j y_j \).  Deriving them from the node locations
    !! in this way accommodates an arbitrary distribution of step sizes.
    real(real64), intent(in), dimension(:) :: t
        !! The interpolation nodes, the first of which is the new point.
    real(real64), intent(out), dimension(:) :: c
        !! An array where the coefficients will be written.

    ! Local Variables
    integer(int32) :: j, m, np
    real(real64) :: p

    ! Process
    np = size(t)
    c(1) = 0.0d0
    do m = 2, np
        c(1) = c(1) + 1.0d0 / (t(1) - t(m))
    end do
    do j = 2, np
        p = 1.0d0
        do m = 2, np
            if (m == j) cycle
            p = p * (t(1) - t(m)) / (t(j) - t(m))
        end do
        c(j) = p / (t(j) - t(1))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine bdf_build_linear_system(this, sys, h, x, y, A, b, args)
    !! Builds the linearized system used to compute the solution for the
    !! implicit step via Newton iteration.
    !!
    !! The corrector residual is
    !! \[
    !! R = M \left( c_0 y + s \right) - f(x,y),
    !! \]
    !! where \(s\) collects the contribution of the stored history, and its
    !! Jacobian is \( c_0 M - J \).  The mass matrix appears only as a factor
    !! and is never inverted, which is what allows a singular \(M\) to be used.
    class(bdf), intent(inout) :: this
        !! The bdf object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.  The step size is already carried by the
        !! coefficients formed by the predictor, so it is not used here.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! The current value(s) of the dependent variable(s).
    real(real64), intent(out), dimension(:,:) :: A
        !! The N-by-N matrix \(A\), where \(A = \frac{\alpha_o}{h_n} M - J\).
    real(real64), intent(out), dimension(:) :: b
        !! The N-element right-hand-side vector for the linear system.
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information in and out
        !! of the differential equation subroutine.

    ! Local Variables
    integer(int32) :: i, n
    logical :: usemass

    ! Initialization
    n = size(y)
    call this%initialize_workspace(n)
    usemass = associated(sys%mass_matrix)

    ! Input Checking
    if (size(A, 1) /= n .or. size(A, 2) /= n) error stop DIFFEQ_MATRIX_SIZE_ERROR
    if (size(b) /= n) error stop DIFFEQ_ARRAY_SIZE_ERROR

    ! Evaluate the model
    call sys%compute_jacobian(x, y, this%m_jac, args)
    call sys%fcn(x, y, this%m_rhs, args)
    this%m_ydot = this%m_c0 * y + this%m_histTerm

    ! Form the system.  Only the derivative of f contributes to the iteration
    ! matrix; the dependence of a state-dependent mass matrix upon y is
    ! neglected, which affects the rate of convergence but not the solution
    ! the iteration converges upon.
    if (usemass) then
        if (.not.this%m_massComputed .or. &
            sys%get_is_mass_matrix_dependent()) &
        then
            call sys%mass_matrix(x, y, this%m_mass, args)
            this%m_massComputed = .true.
        end if
        A = this%m_c0 * this%m_mass - this%m_jac
        b = this%m_rhs - matmul(this%m_mass, this%m_ydot)
    else
        A = -this%m_jac
        do i = 1, n
            A(i,i) = A(i,i) + this%m_c0
        end do
        b = this%m_rhs - this%m_ydot
    end if
end subroutine

! ------------------------------------------------------------------------------
end module