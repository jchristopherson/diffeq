module diffeq_implicit_runge_kutta
    !! Linearly implicit Rosenbrock ODE solvers.
    !!
    !! The solver targets stiff initial-value problems and supports systems
    !! written with a mass matrix,
    !! \[
    !! M(x,y)y' = f(x,y).
    !! \]
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use linalg_qr, only : solve_qr, qr_factor, identity
    implicit none
    private
    public :: rosenbrock
    public :: kennedy_carpenter_4
    public :: kennedy_carpenter_5

    type, abstract, extends(single_step_integrator) :: kennedy_carpenter
        !! Shared implementation for Kennedy--Carpenter ESDIRK methods.
        !!
        !! The stage equations are solved with Newton iteration.  A supplied
        !! mass matrix is incorporated in each stage system as
        !! \(M-h a_{ii}J\), while the embedded tableau supplies the error
        !! estimate used by the inherited adaptive driver.
    contains
        procedure, public :: pre_step_action => kc_pre_step
        procedure, public :: attempt_step => kc_attempt_step
        procedure, public :: post_step_action => kc_post_step
        procedure, public :: interpolate => kc_interpolate
        procedure, public :: get_is_fsal => kc_get_is_fsal
        procedure, public :: get_stage_count => kc_get_stage_count
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_4
        !! Fourth-order ARK4(3)6L[2]SA ESDIRK method.
    contains
        procedure, public :: get_order => kc4_get_order
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_5
        !! Fifth-order ARK5(4)8L[2]SA ESDIRK method.
    contains
        procedure, public :: get_order => kc5_get_order
    end type

    type, extends(single_step_integrator) :: rosenbrock
        !! A fourth-order Rosenbrock integrator for stiff systems.
        !!
        !! Instead of solving a nonlinear equation at every stage, the method
        !! reuses one factored linear system based on the Jacobian
        !! \( J = \partial f / \partial y \).  Without a mass matrix, the
        !! factored matrix is
        !! \[
        !! A = \frac{1}{\gamma h}I - J,
        !! \]
        !! while a user-supplied mass matrix produces
        !! \[
        !! A = \frac{1}{\gamma h}M - J.
        !! \]
        !! The method is appropriate for stiff problems where explicit
        !! Runge--Kutta methods would require prohibitively small steps.
        !! A mass matrix may be constant or state-dependent; state-dependent
        !! matrices are recomputed as needed.
        real(real64), private, allocatable, dimension(:,:) :: jac
            ! The Jacobian matrix.
        real(real64), private, allocatable, dimension(:,:) :: mass
            ! The mass matrix.
        integer(int32), private, allocatable, dimension(:) :: pivot
            ! QR factorization pivot tracking array
        real(real64), private, allocatable, dimension(:) :: tau
            ! QR factorization scalar factors array
        logical, private :: m_massComputed = .false.
            ! True if the mass matrix has been computed; else, false.
        real(real64), private, allocatable, dimension(:) :: dfdx
            ! N-element array of df/dx
        real(real64), private, allocatable, dimension(:,:) :: a
            ! System matrix.
        real(real64), private, allocatable, dimension(:,:) :: qr
            ! QR factored matrix
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
        real(real64), private, allocatable, dimension(:) :: rc4
        logical, private :: m_firstStep = .true.
        logical, private :: m_rejectStep = .false.
        real(real64), private :: m_hOld = 0.0d0
    contains
        procedure, private :: initialize_matrices => rbrk_init_matrices
            !! Allocates internal storage for the system matrices.
        procedure, private :: initialize_interp => rbrk_init_interp
            !! Initializes private storage for the interpolation process.
        procedure, public :: pre_step_action => rbrk_form_matrix
            !! Constructs the system matrix.
        procedure, public :: attempt_step => rbrk_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rbrk_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rbrk_interp
            !! Performs the interpolation.
        procedure, public :: get_order => rbrk_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rbrk_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => rbrk_get_stage_count
            !! Gets the stage count for this integrator.
        procedure, public :: estimate_next_step_size => rbrk_next_step
            !! Estimates the next step size.
    end type

contains
! ******************************************************************************
! ROSENBROCK
! ------------------------------------------------------------------------------
subroutine rbrk_form_matrix(this, prevs, sys, h, x, y, f, args)
    use diffeq_rosenbrock_constants
    !! Constructs the system matrix of the form 
    !! \( A = \frac{1}{\gamma h} M - J \), and then computes it's LU 
    !! factorization.  The LU-factored form of A is stored internally.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
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
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer(int32) :: i, n, flag
    logical :: useMass
    real(real64) :: fac
    real(real64), allocatable, dimension(:,:) :: lu
    
    ! Initialization
    n = size(y)
    useMass = associated(sys%mass_matrix)
    call this%initialize_matrices(n, useMass)

    ! Compute the Jacobian matrix - only need to update if the previous step
    ! was successful
    if (prevs) then
        call sys%compute_jacobian(x, y, this%jac, args = args)

        ! Compute the mass matrix
        if (useMass) then
            if (.not.this%m_massComputed .or. &
                sys%get_is_mass_matrix_dependent()) &
            then
                ! We need to compute the mass matrix
                call sys%mass_matrix(x, y, this%mass, args)
                this%m_massComputed = .true.
            end if
        end if
    end if

    ! Form the system matrix, and then factor it accordingly
    fac = 1.0d0 / (gam * h)
    if (useMass) then
        this%a = fac * this%mass - this%jac
    else
        this%a = -this%jac
        do i = 1, n
            this%a(i,i) = this%a(i,i) + fac
        end do
    end if

    ! Factor the equations
    this%pivot = 0
    call qr_factor(this%a, this%pivot, tau = this%tau, qr = this%qr)

    ! Compute df/dx
    fac = sys%get_finite_difference_step()
    call sys%fcn(x + fac, y, this%dfdx, args)
    this%dfdx = (this%dfdx - f) / h ! forward differencing
end subroutine

! ------------------------------------------------------------------------------
subroutine rbrk_init_matrices(this, n, usemass)
    !! Allocates internal storage for the system matrices.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
    integer(int32), intent(in) :: n
        !! The number of equations being integrated.
    logical, intent(in) :: usemass
        !! True if a mass matrix is used; else, false.

    ! Process
    if (allocated(this%jac)) then
        if (size(this%jac, 1) == n .and. size(this%jac, 2) == n) then
            if (.not.usemass .and. allocated(this%mass)) then
                deallocate(this%mass)
                this%m_massComputed = .false.
            end if
            ! All is good
            return
        else
            deallocate(this%jac)
            if (allocated(this%mass)) deallocate(this%mass)
            deallocate(this%pivot)
            deallocate(this%tau)
            deallocate(this%dfdx)
            deallocate(this%a)
            deallocate(this%qr)
        end if
    end if
    if (usemass) then
        allocate( &
            this%jac(n, n), &
            this%mass(n, n), &
            this%pivot(n), &
            this%tau(n), &
            this%dfdx(n), &
            this%a(n, n), &
            this%qr(n, n) &
        )
    else
        allocate( &
            this%jac(n, n), &
            this%pivot(n), &
            this%tau(n), &
            this%dfdx(n), &
            this%a(n, n), &
            this%qr(n, n) &
        )
    end if

end subroutine

! ------------------------------------------------------------------------------
subroutine rbrk_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    use diffeq_rosenbrock_constants
    !! Attempts an integration step for this integrator.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
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
    integer(int32) :: n
    real(real64), allocatable, dimension(:) :: rhs

    ! Initialization
    n = size(y)
    allocate(rhs(n))

    ! Process
    rhs = f + h * d1 * this%dfdx
    k(:,1) = solve_qr(this%qr, this%tau, this%pivot, rhs)

    yn = y + a21 * k(:,1)
    call sys%fcn(x + c2 * h, yn, fn, args)

    rhs = c21 * k(:,1) / h
    if (allocated(this%mass)) rhs = matmul(this%mass, rhs)
    rhs = fn + h * d2 * this%dfdx + rhs
    k(:,2) = solve_qr(this%qr, this%tau, this%pivot, rhs)

    yn = y + a31 * k(:,1) + a32 * k(:,2)
    call sys%fcn(x + c3 * h, yn, fn, args)

    rhs = (c31 * k(:,1) + c32 * k(:,2)) / h
    if (allocated(this%mass)) rhs = matmul(this%mass, rhs)
    rhs = fn + h * d3 * this%dfdx + rhs
    k(:,3) = solve_qr(this%qr, this%tau, this%pivot, rhs)

    yn = y + a41 * k(:,1) + a42 * k(:,2) + a43 * k(:,3)
    call sys%fcn(x + c4 * h, yn, fn, args)

    rhs = (c41 * k(:,1) + c42 * k(:,2) + c43 * k(:,3)) / h
    if (allocated(this%mass)) rhs = matmul(this%mass, rhs)
    rhs = fn + h * d4 * this%dfdx + rhs
    k(:,4) = solve_qr(this%qr, this%tau, this%pivot, rhs)

    yn = y + a51 * k(:,1) + a52 * k(:,2) + a53 * k(:,3) + a54 * k(:,4)
    call sys%fcn(x + h, yn, fn, args)

    rhs = (c51 * k(:,1) + c52 * k(:,2) + c53 * k(:,3) + &
        c54 * k(:,4)) / h
    if (allocated(this%mass)) rhs = matmul(this%mass, rhs)
    rhs = fn + rhs
    k(:,5) = solve_qr(this%qr, this%tau, this%pivot, rhs)


    yn = yn + k(:,5)
    call sys%fcn(x + h, yn, fn, args) ! updated derivative

    rhs = (c61 * k(:,1) + c62 * k(:,2) + c63 * k(:,3) + &
        c64 * k(:,4) + c65 * k(:,5)) / h
    if (allocated(this%mass)) rhs = matmul(this%mass, rhs)
    yerr = fn + rhs
    yerr = solve_qr(this%qr, this%tau, this%pivot, yerr)

    yn = yn + yerr
end subroutine

! ------------------------------------------------------------------------------
subroutine rbrk_init_interp(this, neqn)
    !! Allocates storage for the interpolation process.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
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
        end if
    end if

    allocate( &
        this%rc1(neqn), &
        this%rc2(neqn), &
        this%rc3(neqn), &
        this%rc4(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rbrk_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    use diffeq_rosenbrock_constants
    !! Sets up the interpolation process.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
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

    ! Local Variables
    integer(int32) :: i, n

    ! Quick Return
    if (.not.dense) return

    ! Initialization
    n = size(y)
    call this%initialize_interp(n)

    ! Process
    do i = 1, n
        this%rc1(i) = y(i)
        this%rc2(i) = yn(i)
        this%rc3(i) = d21 * k(i,1) + d22 * k(i,2) + d23 * k(i,3) + &
            d24 * k(i,4) + d25 * k(i,5)
        this%rc4(i) = d31 * k(i,1) + d32 * k(i,2) + d33 * k(i,3) + &
            d34 * k(i,4) + d35 * k(i,5)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine rbrk_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(rosenbrock), intent(in) :: this
        !! The rosenbrock object.
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
    y = this%rc1 * s1 + s * (this%rc2 + s1 * (this%rc3 + s * this%rc4))
end subroutine

! ------------------------------------------------------------------------------
pure function rbrk_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(rosenbrock), intent(in) :: this
        !! The rosenbrock object.
    integer(int32) :: rst
        !! The order.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function rbrk_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(rosenbrock), intent(in) :: this
        !! The rosenbrock object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function rbrk_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(rosenbrock), intent(in) :: this
        !! The rosenbrock object.
    integer(int32) :: rst
        !! The stage count.
    rst = 6
end function

! ------------------------------------------------------------------------------
function rbrk_next_step(this, e, eold, h, x) result(rst)
    !! Estimates the next step size based upon the current and previous error
    !! estimates.
    class(rosenbrock), intent(inout) :: this
        !! The rosenbrock object.
    real(real64), intent(in) :: e
        !! The norm of the current scaled error estimate.
    real(real64), intent(inout) :: eold
        !! The norm of the previous step's scaled error estimate.  On output,
        !! this variable is updated.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current independent variable value.
    real(real64) :: rst
        !! The new step size estimate.

    ! Parameters
    real(real64), parameter :: fac1 = 5.0d0
    real(real64), parameter :: fac2 = 1.0d0 / 6.0d0

    ! Local Variables
    real(real64) :: safe, fac, facpred

    ! Initialization
    safe = this%get_step_size_factor()
    fac = max(fac2, min(fac1, e**(0.25d0) / safe))
    rst = h / fac

    ! Process
    if (e <= 1.0d0) then
        ! The step met error tolerances and is acceptable
        if (.not.this%m_firstStep) then
            facpred = (this%m_hOld / h) * (e * e / eold)**(0.25d00) / safe
            facpred = max(fac2, min(fac1, facpred))
            fac = max(fac, facpred)
            rst = h / fac
        end if
        this%m_firstStep = .false.
        this%m_hOld = h
        eold = max(1.0d-2, e)
        if (this%m_rejectStep) then
            ! Don't let the step size increase if the last step was rejected
            if (abs(h) >= 0.0d0) then
                rst = min(abs(rst), abs(h))
            else
                rst = max(abs(rst), abs(h))
            end if
            rst = sign(rst, h)
        end if
        this%m_rejectStep = .false.
    else
        ! The step is rejected, reduce the step size - already computed
        this%m_rejectStep = .true.
    end if
end function

! ******************************************************************************
! KENNEDY-CARPENTER INTEGRATORS
! ------------------------------------------------------------------------------
subroutine kc_pre_step(this, prevs, sys, h, x, y, f, args)
    class(kennedy_carpenter), intent(inout) :: this
    logical, intent(in) :: prevs
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x, y(:), f(:)
    class(*), intent(inout), optional :: args
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_post_step(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    class(kennedy_carpenter), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    logical, intent(in) :: dense
    real(real64), intent(in) :: x, xn, y(:), yn(:), f(:), fn(:)
    real(real64), intent(inout) :: k(:,:)
    class(*), intent(inout), optional :: args
end subroutine

! ------------------------------------------------------------------------------
pure function kc_get_is_fsal(this) result(rst)
    class(kennedy_carpenter), intent(in) :: this
    logical :: rst
    rst = .false.
end function

! ------------------------------------------------------------------------------
pure function kc_get_stage_count(this) result(rst)
    class(kennedy_carpenter), intent(in) :: this
    integer(int32) :: rst
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
    class(kennedy_carpenter_4), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function kc5_get_order(this) result(rst)
    class(kennedy_carpenter_5), intent(in) :: this
    integer(int32) :: rst
    rst = 5
end function

! ------------------------------------------------------------------------------
subroutine kc_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    !! Attempts one Kennedy--Carpenter ESDIRK step.
    !!
    !! The stage state is solved from
    !! \[
    !! M(Y_i-Y_i^*)-h a_{ii}f(x+c_i h,Y_i)=0,
    !! \]
    !! using Newton iteration.  The high- and embedded-order solutions are
    !! \(y_{n+1}=y_n+h\sum b_i k_i\) and
    !! \(\widehat y_{n+1}=y_n+h\sum d_i k_i\), respectively.
    class(kennedy_carpenter), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x, y(:), f(:)
    real(real64), intent(out) :: yn(:), fn(:), yerr(:)
    real(real64), intent(out) :: k(:,:)
    class(*), intent(inout), optional :: args
    real(real64) :: a(8,8), b(8), d(8), c(8)
    real(real64) :: base(size(y)), state(size(y)), rhs(size(y)), initial_f(size(y))
    real(real64) :: jac(size(y),size(y)), mass(size(y),size(y))
    real(real64) :: delta(size(y))
    integer(int32) :: stages, i, j, iteration
    logical :: use_mass

    call kc_table(this, a, b, d, c, stages)
    k = 0.0d0
    initial_f = f
    if (associated(sys%mass_matrix)) then
        call sys%mass_matrix(x, y, mass, args)
        call solve_kc_system(mass, initial_f, k(:,1))
    else
        k(:,1) = initial_f
    end if
    do i = 2, stages
        base = y
        do j = 1, i - 1
            base = base + h * a(i,j) * k(:,j)
        end do
        call sys%fcn(x + c(i)*h, base, rhs, args)
        if (associated(sys%mass_matrix)) then
            call sys%mass_matrix(x + c(i)*h, base, mass, args)
            call solve_kc_system(mass, rhs, initial_f)
            state = base + h*a(i,i)*initial_f
        else
            state = base + h*a(i,i)*rhs
        end if
        do iteration = 1, 12
            call sys%fcn(x + c(i)*h, state, rhs, args)
            use_mass = associated(sys%mass_matrix)
            if (use_mass) then
                call sys%mass_matrix(x + c(i)*h, state, mass, args)
                call sys%compute_jacobian(x + c(i)*h, state, jac, args)
                rhs = matmul(mass, (state - base)) - h*a(i,i)*rhs
                    call solve_kc_system(mass - h*a(i,i)*jac, -rhs, delta)
            else
                call sys%compute_jacobian(x + c(i)*h, state, jac, args)
                rhs = state - base - h*a(i,i)*rhs
                    call solve_kc_system(identity(size(y)) - h*a(i,i)*jac, -rhs, delta)
            end if
            if (norm2(rhs) <= 1.0d-12*max(1.0d0,norm2(state))) exit
            state = state + delta
        end do
        if (iteration > 12) error stop DIFFEQ_CONVERGENCE_ERROR
        call sys%fcn(x + c(i)*h, state, rhs, args)
        if (associated(sys%mass_matrix)) then
            call sys%mass_matrix(x + c(i)*h, state, mass, args)
            call solve_kc_system(mass, rhs, k(:,i))
        else
            k(:,i) = rhs
        end if
    end do
    yn = y
    rhs = 0.0d0
    do i = 1, stages
        yn = yn + h*b(i)*k(:,i)
        rhs = rhs + h*d(i)*k(:,i)
    end do
    yerr = yn - (y + rhs)
    call sys%fcn(x + h, yn, fn, args)
    if (associated(sys%mass_matrix)) then
        call sys%mass_matrix(x + h, yn, mass, args)
        call solve_kc_system(mass, fn, initial_f)
        fn = initial_f
    end if
end subroutine

subroutine kc_table(this, a, b, d, c, stages)
    class(kennedy_carpenter), intent(in) :: this
    real(real64), intent(out) :: a(8,8), b(8), d(8), c(8)
    integer(int32), intent(out) :: stages
    a = 0.0d0; b = 0.0d0; d = 0.0d0; c = 0.0d0
    select type (this)
    type is (kennedy_carpenter_4)
        stages = 6
        a(2,1)=1.0d0/4.0d0; a(2,2)=1.0d0/4.0d0
        a(3,1)=8611.0d0/62500.0d0; a(3,2)=-1743.0d0/31250.0d0; a(3,3)=1.0d0/4.0d0
        a(4,1)=5012029.0d0/34652500.0d0; a(4,2)=-654441.0d0/2922500.0d0
        a(4,3)=174375.0d0/388108.0d0; a(4,4)=1.0d0/4.0d0
        a(5,1)=15267082809.0d0/155376265600.0d0; a(5,2)=-71443401.0d0/120774400.0d0
        a(5,3)=730878875.0d0/902184768.0d0; a(5,4)=2285395.0d0/8070912.0d0; a(5,5)=1.0d0/4.0d0
        a(6,1)=82889.0d0/524892.0d0; a(6,3)=15625.0d0/83664.0d0
        a(6,4)=69875.0d0/102672.0d0; a(6,5)=-2260.0d0/8211.0d0; a(6,6)=1.0d0/4.0d0
        b(1)=a(6,1); b(3)=a(6,3); b(4)=a(6,4); b(5)=a(6,5); b(6)=a(6,6)
        d(1)=4586570599.0d0/29645900160.0d0; d(3)=178811875.0d0/945068544.0d0
        d(4)=814220225.0d0/1159782912.0d0; d(5)=-3700637.0d0/11593932.0d0; d(6)=61727.0d0/225920.0d0
        c(2)=0.5d0; c(3)=83.0d0/250.0d0; c(4)=31.0d0/50.0d0; c(5)=17.0d0/20.0d0; c(6)=1.0d0
    type is (kennedy_carpenter_5)
        stages = 8
        a(2,1)=41.0d0/200.0d0; a(2,2)=41.0d0/200.0d0
        a(3,1)=41.0d0/400.0d0; a(3,2)=-567603406766.0d0/11931857230679.0d0; a(3,3)=41.0d0/200.0d0
        a(4,1)=683785636431.0d0/9252920307686.0d0; a(4,3)=-110385047103.0d0/1367015193373.0d0; a(4,4)=41.0d0/200.0d0
        a(5,1)=3016520224154.0d0/10081342136671.0d0; a(5,3)=30586259806659.0d0/12414158314087.0d0
        a(5,4)=-22760509404356.0d0/11113319521817.0d0; a(5,5)=41.0d0/200.0d0
        a(6,1)=218866479029.0d0/1489978393911.0d0; a(6,3)=638256894668.0d0/5436446318841.0d0
        a(6,4)=-1179710474555.0d0/5321154724896.0d0; a(6,5)=-60928119172.0d0/8023461067671.0d0; a(6,6)=41.0d0/200.0d0
        a(7,1)=1020004230633.0d0/5715676835656.0d0; a(7,3)=25762820946817.0d0/25263940353407.0d0
        a(7,4)=-2161375909145.0d0/9755907335909.0d0; a(7,5)=-211217309593.0d0/5846859502534.0d0
        a(7,6)=-4269925059573.0d0/7827059040749.0d0; a(7,7)=41.0d0/200.0d0
        a(8,1)=-872700587467.0d0/9133579230613.0d0; a(8,4)=22348218063261.0d0/9555858737531.0d0
        a(8,5)=-1143369518992.0d0/8141816002931.0d0; a(8,6)=-39379526789629.0d0/19018526304540.0d0
        a(8,7)=32727382324388.0d0/42900044865799.0d0; a(8,8)=41.0d0/200.0d0
        b(1)=a(8,1); b(4)=a(8,4); b(5)=a(8,5); b(6)=a(8,6); b(7)=a(8,7); b(8)=a(8,8)
        d(1)=-975461918565.0d0/9796059967033.0d0; d(4)=78070527104295.0d0/32432590147079.0d0
        d(5)=-548382580838.0d0/3424219808633.0d0; d(6)=-33438840321285.0d0/15594753105479.0d0
        d(7)=3629800801594.0d0/4656183773603.0d0; d(8)=4035322873751.0d0/18575991585200.0d0
        c(2)=41.0d0/100.0d0; c(3)=2935347310677.0d0/11292855782101.0d0; c(4)=1426016391358.0d0/7196633302097.0d0
        c(5)=0.92d0; c(6)=0.24d0; c(7)=0.6d0; c(8)=1.0d0
    end select
end subroutine

subroutine solve_kc_system(a, b, x)
    real(real64), intent(in) :: a(:,:), b(:)
    real(real64), intent(out) :: x(:)

    integer(int32) :: n
    integer(int32), allocatable, dimension(:) :: pivot
    real(real64), allocatable, dimension(:) :: tau
    real(real64), allocatable, dimension(:,:) :: qr

    n = size(a, 2)
    allocate(pivot(n), source = 0)
    call qr_factor(a, pivot, tau = tau, qr = qr)
    x = solve_qr(qr, tau, pivot, b)
end subroutine

subroutine kc_interpolate(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    class(kennedy_carpenter), intent(in) :: this
    real(real64), intent(in) :: x, xn, yn(:), fn(:), xn1, yn1(:), fn1(:)
    real(real64), intent(out) :: y(:)
    real(real64) :: s, h
    h=xn1-xn; s=(x-xn)/h
    y=(2.0d0*s**3-3.0d0*s**2+1.0d0)*yn+(s**3-2.0d0*s**2+s)*h*fn &
        +(-2.0d0*s**3+3.0d0*s**2)*yn1+(s**3-s**2)*h*fn1
end subroutine

! ------------------------------------------------------------------------------
end module