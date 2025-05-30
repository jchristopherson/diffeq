module diffeq_implicit_runge_kutta
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use linalg, only : lu_factor, solve_lu
    implicit none
    private
    public :: rosenbrock

    type, extends(single_step_integrator) :: rosenbrock
        !! Defines a 4th order Rosenbrock integrator.
        !!
        !! Remarks:
        !! 1. This integrator is suitable for systems of stiff equations 
        !!  with modest accuracy requirements.
        !! 2. This integrator is capable of dealing with systems that utilize a 
        !!  mass matrix.
        real(real64), private, allocatable, dimension(:,:) :: jac
            ! The Jacobian matrix.
        real(real64), private, allocatable, dimension(:,:) :: mass
            ! The mass matrix.
        integer(int32), private, allocatable, dimension(:) :: pivot
            ! LU factorization pivot tracking array
        logical, private :: m_massComputed = .false.
            ! True if the mass matrix has been computed; else, false.
        real(real64), private, allocatable, dimension(:) :: dfdx
            ! N-element array of df/dx
        real(real64), private, allocatable, dimension(:,:) :: a
            ! System matrix.
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
subroutine rbrk_form_matrix(this, prevs, sys, h, x, y, f, args, err)
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.

    ! Local Variables
    integer(int32) :: i, n
    logical :: useMass
    real(real64) :: fac
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(y)
    useMass = associated(sys%mass_matrix)
    call this%initialize_matrices(n, useMass)

    ! Compute the Jacobian matrix - only need to update if the previous step
    ! was successful
    if (prevs) then
        call sys%compute_jacobian(x, y, this%jac, args = args, err = errmgr)
        if (errmgr%has_error_occurred()) return

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
    call lu_factor(this%a, this%pivot, errmgr)
    if (errmgr%has_error_occurred()) return

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
            ! All is good
            return
        else
            deallocate(this%jac)
            if (allocated(this%mass)) deallocate(this%mass)
            deallocate(this%pivot)
            deallocate(this%dfdx)
            deallocate(this%a)
        end if
    end if
    if (usemass) then
        allocate( &
            this%jac(n, n), &
            this%mass(n, n), &
            this%pivot(n), &
            this%dfdx(n), &
            this%a(n, n) &
        )
    else
        allocate( &
            this%jac(n, n), &
            this%pivot(n), &
            this%dfdx(n), &
            this%a(n, n) &
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

    ! Initialization
    n = size(y)

    ! Process
    k(:,1) = f + h * d1 * this%dfdx
    call solve_lu(this%a, this%pivot, k(:,1))

    yn = y + a21 * k(:,1)
    call sys%fcn(x + c2 * h, yn, fn, args)

    k(:,2) = fn + h * d2 * this%dfdx + c21 * k(:,1) / h
    call solve_lu(this%a, this%pivot, k(:,2))

    yn = y + a31 * k(:,1) + a32 * k(:,2)
    call sys%fcn(x + c3 * h, yn, fn, args)

    k(:,3) = fn + h * d3 * this%dfdx + (c31 * k(:,1) + c32 * k(:,2)) / h
    call solve_lu(this%a, this%pivot, k(:,3))

    yn = y + a41 * k(:,1) + a42 * k(:,2) + a43 * k(:,3)
    call sys%fcn(x + c4 * h, yn, fn, args)

    k(:,4) = fn + h * d4 * this%dfdx + (c41 * k(:,1) + c42 * k(:,2) + &
        c43 * k(:,3)) / h
    call solve_lu(this%a, this%pivot, k(:,4))

    yn = y + a51 * k(:,1) + a52 * k(:,2) + a53 * k(:,3) + a54 * k(:,4)
    call sys%fcn(x + h, yn, fn, args)

    k(:,5) = fn + (c51 * k(:,1) + c52 * k(:,2) + c53 * k(:,3) + &
        c54 * k(:,4)) / h
    call solve_lu(this%a, this%pivot, k(:,5))

    yn = yn + k(:,5)
    call sys%fcn(x + h, yn, fn, args) ! updated derivative

    yerr = fn + (c61 * k(:,1) + c62 * k(:,2) + c63 * k(:,3) + &
        c64 * k(:,4) + c65 * k(:,5)) / h
    call solve_lu(this%a, this%pivot, yerr)

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
function rbrk_next_step(this, e, eold, h, x, err) result(rst)
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
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.
        !!
        !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size
        !!      becomes too small in magnitude.
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

! ------------------------------------------------------------------------------
end module