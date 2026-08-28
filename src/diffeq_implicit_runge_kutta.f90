module diffeq_implicit_runge_kutta
    !! Implicit Runge-Kutta ODE solvers.
    !!
    !! The solvers target stiff initial-value problems and support systems
    !! written with a mass matrix,
    !! \[
    !! M(x,y)y' = f(x,y).
    !! \]
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use linalg, only : solve_qr, qr_factor, identity, svd
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
        !! estimate used by the inherited adaptive driver.  These solvers use
        !! the inherited configurable PI step-size controller.
    contains
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
        !! matrices are recomputed as needed.  Step sizes are selected by the
        !! Rosenbrock-specific controller rather than the inherited PI
        !! controller.
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
            ! Interpolation coefficient array.
        real(real64), private, allocatable, dimension(:) :: rc2
            ! Interpolation coefficient array.
        real(real64), private, allocatable, dimension(:) :: rc3
            ! Interpolation coefficient array.
        real(real64), private, allocatable, dimension(:) :: rc4
            ! Interpolation coefficient array.
        logical, private :: m_firstStep = .true.
            ! True if the integrator is on its first step; else, false.
        logical, private :: m_rejectStep = .false.
            ! True if the previous step was rejected; else, false.
        real(real64), private :: m_hOld = 0.0d0
            ! The previous step size.
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
    !! \( A = \frac{1}{\gamma h} M - J \), and then computes its QR
    !! factorization.  The QR-factored form of A is stored internally.
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
    !! Placeholder routine for any pre-step actions.
    !!
    !! The Jacobian and mass matrices are state-dependent within each stage,
    !! so they are formed by the stage iteration rather than here.
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
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_post_step(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    !! Placeholder routine for any post-step actions.
    !!
    !! The interpolation for these integrators requires only the solution and
    !! derivative values at each end of the step, so no additional storage is
    !! required here.
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

    ! Process
    return
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
    integer(int32) :: rst
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
    integer(int32) :: rst
        !! The order.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function kc5_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(kennedy_carpenter_5), intent(in) :: this
        !! The kennedy_carpenter_5 object.
    integer(int32) :: rst
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

    ! Parameters
    integer(int32), parameter :: maxiter = 12
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    logical :: usemass
    integer(int32) :: i, j, n, stages, iteration
    real(real64) :: a(8,8), b(8), d(8), c(8)
    real(real64), allocatable, dimension(:) :: base, state, rhs, deriv, &
        residual, delta, embedded
    real(real64), allocatable, dimension(:,:) :: jac, mass

    ! Initialization
    n = size(y)
    usemass = associated(sys%mass_matrix)
    call kc_table(this, a, b, d, c, stages)
    allocate( &
        base(n), &
        state(n), &
        rhs(n), &
        deriv(n), &
        residual(n), &
        delta(n), &
        embedded(n), &
        jac(n, n), &
        mass(n, n) &
    )

    ! Process
    ! The first stage of an ESDIRK method is explicit
    k = 0.0d0
    if (usemass) then
        call sys%mass_matrix(x, y, mass, args)
        call kc_consistent_derivative(sys, x, y, f, mass, k(:,1), args)
    else
        k(:,1) = f
    end if

    ! Each remaining stage is diagonally implicit
    do i = 2, stages
        ! Accumulate the contributions from the previously computed stages
        base = y
        do j = 1, i - 1
            base = base + h * a(i,j) * k(:,j)
        end do

        if (usemass) then
            call kc_mass_stage(sys, x + c(i) * h, base, h * a(i,i), &
                state, k(:,i), args)
            cycle
        end if

        ! Use an explicit Euler prediction to start the Newton iteration
        call sys%fcn(x + c(i) * h, base, rhs, args)
        state = base + h * a(i,i) * rhs

        ! Solve the stage equation
        do iteration = 1, maxiter
            call sys%fcn(x + c(i) * h, state, rhs, args)
            if (usemass) then
                call sys%mass_matrix(x + c(i) * h, state, mass, args)
                call sys%compute_jacobian(x + c(i) * h, state, jac, args)
                residual = matmul(mass, state - base) - h * a(i,i) * rhs
                delta = solve_kc_system(mass - h * a(i,i) * jac, -residual)
            else
                call sys%compute_jacobian(x + c(i) * h, state, jac, args)
                residual = state - base - h * a(i,i) * rhs
                delta = solve_kc_system(identity(n) - h * a(i,i) * jac, &
                    -residual)
            end if
            if (norm2(residual) <= tol * max(1.0d0, norm2(state))) exit
            state = state + delta
        end do
        if (iteration > maxiter) error stop DIFFEQ_CONVERGENCE_ERROR

        ! Store the derivative for this stage
        call sys%fcn(x + c(i) * h, state, rhs, args)
        if (usemass) then
            call sys%mass_matrix(x + c(i) * h, state, mass, args)
            k(:,i) = solve_kc_system(mass, rhs)
        else
            k(:,i) = rhs
        end if
    end do

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
subroutine kc_consistent_derivative(sys, x, y, f, mass, derivative, args)
    !! Computes a consistent initial derivative for an index-1 DAE.
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, y(:), f(:), mass(:,:)
    real(real64), intent(out) :: derivative(:)
    class(*), intent(inout), optional :: args
    real(real64), allocatable :: singular_values(:), left_vectors(:,:), &
        jac(:,:), fplus(:), augmented(:,:), rhs(:)
    real(real64) :: fdstep, rank_tolerance, scale
    integer(int32) :: n, rank, nullity, i

    n = size(y)
    call svd(mass, singular_values, left_vectors)
    scale = max(1.0d0, singular_values(1))
    rank_tolerance = 100.0d0 * epsilon(1.0d0) * scale
    rank = count(singular_values > rank_tolerance)

    if (rank == n) then
        derivative = solve_kc_system(mass, f)
        return
    end if

    nullity = n - rank
    allocate(jac(n,n), fplus(n), augmented(n,n), rhs(n))
    call sys%compute_jacobian(x, y, jac, args)
    fdstep = sys%get_finite_difference_step()
    call sys%fcn(x + fdstep, y, fplus, args)
    fplus = (fplus - f) / fdstep

    augmented = 0.0d0
    rhs = 0.0d0
    augmented(1:rank,:) = matmul(transpose(left_vectors(:,1:rank)), mass)
    rhs(1:rank) = matmul(transpose(left_vectors(:,1:rank)), f)
    do i = 1, nullity
        augmented(rank+i,:) = matmul(left_vectors(:,rank+i), jac)
        rhs(rank+i) = -dot_product(left_vectors(:,rank+i), fplus)
    end do

    derivative = solve_kc_system(augmented, rhs)
end subroutine

subroutine kc_mass_stage(sys, x, base, diagonal_step, state, derivative, args)
    !! Solves one Kennedy--Carpenter mass-matrix stage without forming M^{-1}f.
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, base(:), diagonal_step
    real(real64), intent(out) :: state(:), derivative(:)
    class(*), intent(inout), optional :: args
    real(real64) :: f(size(base)), residual(size(base))
    real(real64) :: jac(size(base),size(base)), mass(size(base),size(base))
    real(real64) :: mass_perturbed(size(base),size(base))
    real(real64) :: system(size(base),size(base)), delta(size(base))
    real(real64) :: fdstep
    integer(int32) :: i, iteration

    state = base
    do iteration = 1, 12
        derivative = (state - base) / diagonal_step
        call sys%fcn(x, state, f, args)
        call sys%mass_matrix(x, state, mass, args)
        residual = matmul(mass, state - base) - diagonal_step * f
        if (norm2(residual) <= 1.0d-12 * max(1.0d0, norm2(state))) exit
        call sys%compute_jacobian(x, state, jac, args)
        system = mass - diagonal_step * jac
        fdstep = sys%get_finite_difference_step()
        do i = 1, size(base)
            state(i) = state(i) + fdstep
            call sys%mass_matrix(x, state, mass_perturbed, args)
            state(i) = state(i) - fdstep
            system(:,i) = system(:,i) + &
                matmul((mass_perturbed - mass) / fdstep, state - base)
        end do
        delta = solve_kc_system(system, -residual)
        state = state + delta
    end do
    if (iteration > 12) error stop DIFFEQ_CONVERGENCE_ERROR
    derivative = (state - base) / diagonal_step
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
    integer(int32), intent(out) :: stages
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

    ! Local Variables
    integer(int32) :: n, rank
    integer(int32), allocatable, dimension(:) :: pivot
    real(real64), allocatable, dimension(:) :: tau
    real(real64), allocatable, dimension(:,:) :: qr
    real(real64), allocatable, dimension(:) :: singular_values
    real(real64) :: rank_tolerance

    ! Initialization
    n = size(a, 2)
    allocate(pivot(n), source = 0)

    call svd(a, singular_values)
    rank_tolerance = 100.0d0 * epsilon(1.0d0) * &
        max(1.0d0, singular_values(1))
    rank = count(singular_values > rank_tolerance)
    if (rank /= n) error stop DIFFEQ_SINGULAR_MATRIX_ERROR

    ! Process
    call qr_factor(a, pivot, tau = tau, qr = qr)
    x = solve_qr(qr, tau, pivot, b)
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
        (s**3 - 2.0d0 * s**2 + s) * h * fn + &
        (-2.0d0 * s**3 + 3.0d0 * s**2) * yn1 + &
        (s**3 - s**2) * h * fn1
end subroutine

! ------------------------------------------------------------------------------
end module