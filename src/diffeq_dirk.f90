module diffeq_dirk
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use ferror
    implicit none
    private
    public :: dirk_get_matrix_value
    public :: dirk_get_array_value
    public :: dirk_eval_error
    public :: diagonally_implicit_integrator
    public :: implicit_runge_kutta_4
    
    type, abstract, extends(single_step_integrator) :: &
        diagonally_implicit_integrator
        !! Defines a diagonally implicit, single-step integrator.
        logical, private :: m_matrixCurrent = .false.
            ! Is the iteration matrix current
        real(real64), private :: m_tol = 1.0d-8
            ! The Newton iteration convergence tolerance.
        integer(int32), private :: m_maxIter = 10
            ! The maximum number of Newton iterations allowed.
        real(real64), private, allocatable, dimension(:) :: xi
            ! N-element Newton solver workspace array.
        real(real64), private, allocatable, dimension(:) :: dy
            ! N-element Newton solver change in solution array.
    contains
        procedure, private :: initialize_solver => dii_init_solver_storage
            !! Allocates internal resources used by the Newton solver.
        procedure(dirk_get_matrix_value), public, pass, deferred :: &
            get_model_coefficient
            !! Gets the requested model coefficient from the Butcher tableau.
        procedure(dirk_get_array_value), public, pass, deferred :: get_weight
            !! Gets the requested weighting factor from the Butcher tableau.
        procedure(dirk_get_array_value), public, pass, deferred :: get_node
            !! Gets the requested node value from the Butcher tableau.
        procedure(dirk_eval_error), public, pass, deferred :: compute_errors
            !! Computes the error for each equation.
        procedure, public :: get_is_matrix_current => dii_get_mtx_current
            !! Gets a value determining if the iteration matrix is current.
        procedure, public :: set_is_matrix_current => dii_set_mtx_current
            !! Sets a value determining if the iteration matrix is current.
        procedure, public :: form_matrix => dii_form_factored_matrix
            !! Forms the iteration matrix and computes its LU factorization.
        procedure, public :: get_newton_tolerance => dii_get_newton_tol
            !! Gets the convergence tolerance for the Newton solver.
        procedure, public :: set_newton_tolerance => dii_set_newton_tol
            !! Sets the convergence tolerance for the Newton solver.
        procedure, public :: get_newton_iteration_limit => &
            dii_get_newton_iter_limit
            !! Gets the maximum number of Newton iterations allowed.
        procedure, public :: set_newton_iteration_limit => &
            dii_set_newton_iter_limit
            !! Sets the maximum number of Newton iterations allowed.
        procedure, public :: newton_iteration => dii_solve_newton
            !! Performs the Newton iteration for the i-th stage.
        procedure, public :: solve => dii_ode_solver
            !! Solves the supplied system of ODE's.
        procedure, public :: attempt_step => dii_attempt_step
            !! This is a pass-thru routine that has no function for DIRK 
            !! integrators.
    end type

    interface
        pure function dirk_get_matrix_value(this, i, j) result(rst)
            !! Retrieves a matrix value from the integrator.
            use iso_fortran_env
            import diagonally_implicit_integrator
            class(diagonally_implicit_integrator), intent(in) :: this
                !! The diagonally_implicit_integrator object.
            integer(int32), intent(in) :: i
                !! The row index.
            integer(int32), intent(in) :: j
                !! The column index.
            real(real64) :: rst
                !! The value.
        end function

        pure function dirk_get_array_value(this, i) result(rst) 
            !! Retrieves an array value from the integrator.
            use iso_fortran_env
            import diagonally_implicit_integrator
            class(diagonally_implicit_integrator), intent(in) :: this
                !! The diagonally_implicit_integrator object.
            integer(int32), intent(in) :: i
                !! The array index.
            real(real64) :: rst
                !! The value.
        end function

        subroutine dirk_eval_error(this, sys, x, y, f, xn, yn, fn, k, resid, &
            niter, yerr)
            !! Computes the error for each equation.
            use iso_fortran_env
            import diagonally_implicit_integrator
            import ode_container
            class(diagonally_implicit_integrator), intent(in) :: this
                !! The diagonally_implicit_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the solution at x.
            real(real64), intent(in), dimension(:) :: f
                !! An N-element array containing the derivatives at x.
            real(real64), intent(in) :: xn
                !! The next value of the independent variable.
            real(real64), intent(in), dimension(:) :: yn
                !! An N-element array containing the solution estimate at xn.
            real(real64), intent(in), dimension(:) :: fn
                !! An N-element array containing the derivatives at xn.
            real(real64), intent(in), dimension(:,:) :: k
                !! An N-by-NSTAGES matrix containing the derivatives at each
                !! stage.
            real(real64), intent(in) :: resid
                !! The norm of the residual from the Newton iteration process.
            integer(int32), intent(in) :: niter
                !! The number of iterations taken by the Newton process.
            real(real64), intent(out), dimension(:) :: yerr
                !! An N-element array where this routine will write the error
                !! estimates for each equation.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, extends(diagonally_implicit_integrator) :: implicit_runge_kutta_4
        !! A 6-stage, 4th order implicit Runge-Kutta integrator.
        real(real64), private :: m_a(6, 6)
            ! Model coefficients - Butcher's tableau
        real(real64), private :: m_b(6)
            ! Weighting factors
        real(real64), private :: m_c(6)
            ! Nodes
        real(real64), private :: m_dc(4, 6)
            ! Interpolation coefficients
        real(real64), private :: m_e(6)
            ! Error coefficients
        logical, private :: m_populated = .false.
            ! Are the model coefficient matrices populated?
        real(real64), allocatable, dimension(:,:) :: m_k
            ! Derivative values at each stage (N-by-NSTAGE).
    contains
        procedure, public :: pre_step_action => irk4_pre_step
            !! Performs any pre-step actions.
        procedure, public :: get_order => irk4_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => irk4_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => irk4_get_stage_count
            !! Gets the stage count for this integrator.
        procedure, public :: get_model_coefficient => irk4_get_model_coeff
            !! Gets the requested model coefficient from the Butcher tableau.
        procedure, public :: get_weight => irk4_get_weight
            !! Gets the requested weighting factor from the Butcher tableau.
        procedure, public :: get_node => irk4_get_node
            !! Gets the requested weighting factor from the Butcher tableau.
        procedure, public :: compute_errors => irk4_eval_error
            !! Computes the error for each equation.
        procedure, public :: post_step_action => irk4_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => irk4_interp
            !! Performs the interpolation.
    end type

contains

! ******************************************************************************
! DIAGONALLY IMPLICIT INTEGRATOR
! ------------------------------------------------------------------------------
pure function dii_get_mtx_current(this) result(rst)
    !! Gets a value determining if the iteration matrix is current.
    class(diagonally_implicit_integrator), intent(in) :: this
        !! The diagonally_implicit_integrator object.
    logical :: rst
        !! True if the matrix is current; else, false.
    rst = this%m_matrixCurrent
end function

! --------------------
subroutine dii_set_mtx_current(this, x)
    !! Sets a value determining if the iteration matrix is current.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    logical, intent(in) :: x
        !! True if the matrix is current; else, false.
    this%m_matrixCurrent = x
end subroutine

! ------------------------------------------------------------------------------
subroutine dii_form_factored_matrix(this, sys, alpha, h, x, y, jac, a, pvt, err)
    use linalg, only : lu_factor
    !! Forms the iteration matrix of the form \( A = M - h \alpha J \)
    !! or \( A = I - h \alpha J \), and then computes its LU factorization.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the equations to integrate.
    real(real64), intent(in) :: alpha
        !! The scalar multiplier to the Jacobian matrix.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution values at x.
    real(real64), intent(out), dimension(:,:) :: jac
        !! An N-by-N matrix containing the Jacobian.
    real(real64), intent(out), dimension(:,:) :: a
        !! The N-by-N LU-factored iteration matrix.
    integer(int32), intent(out), dimension(:) :: pvt
        !! The N-element pivot tracking array from the LU factorization.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to 
        !! provide error handling.  Possible errors and warning messages
        !! that may be encountered are as follows.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    integer(int32) :: i, n
    logical :: useMass
    real(real64) :: fac
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(y)
    fac = h * alpha

    ! Compute the Jacobian
    call sys%compute_jacobian(x, y, jac, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Form the matrix
    useMass = associated(sys%mass_matrix)
    if (useMass) then
        call sys%mass_matrix(x, y, a)
        a = a - fac * jac
    else
        a = -fac * jac
        do i = 1, n
            a(i,i) = a(i,i) + 1.0d0
        end do
    end if

    ! Compute the LU factorization of the matrix
    call lu_factor(a, pvt, errmgr)
    if (errmgr%has_error_occurred()) return

    ! If we're here, everything is updated
    call this%set_is_matrix_current(.true.)
end subroutine

! ------------------------------------------------------------------------------
subroutine dii_init_solver_storage(this, n)
    !! Allocates internal resources used by the Newton solver.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    integer(int32), intent(in) :: n
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%xi)) then
        if (size(this%xi) == n) then
            ! All is good
            return
        else
            deallocate(this%xi)
            deallocate(this%dy)
        end if
    end if
    allocate( &
        this%xi(n), &
        this%dy(n) &
    )
end subroutine

! ------------------------------------------------------------------------------
pure function dii_get_newton_tol(this) result(rst)
    !! Gets the convergence tolerance for the Newton solver.
    class(diagonally_implicit_integrator), intent(in) :: this
        !! The diagonally_implicit_integrator object.
    real(real64) :: rst
        !! The tolerance.
    rst = this%m_tol
end function

! --------------------
subroutine dii_set_newton_tol(this, x)
    !! Gets the convergence tolerance for the Newton solver.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance.
    this%m_tol = x
end subroutine

! ------------------------------------------------------------------------------
pure function dii_get_newton_iter_limit(this) result(rst)
    !! Gets the maximum number of Newton iterations allowed.
    class(diagonally_implicit_integrator), intent(in) :: this
        !! The diagonally_implicit_integrator object.
    integer(int32) :: rst
        !! The limit.
    rst = this%m_maxIter
end function

! --------------------
subroutine dii_set_newton_iter_limit(this, x)
    !! Sets the maximum number of Newton iterations allowed.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    integer(int32), intent(in) :: x
        !! The limit.
    this%m_maxIter = x
end subroutine

! ------------------------------------------------------------------------------
subroutine dii_solve_newton(this, sys, a, pvt, i, h, x, y, yn, fn, niter, &
    accept, resid)
    use linalg, only : solve_lu
    !! Performs the Newton iteration for the i-th stage.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the equations to integrate.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N LU-factored iteration matrix.
    integer(int32), intent(in), dimension(:) :: pvt
        !! The N-element pivot tracking array from the LU factorization.
    integer(int32), intent(in) :: i
        !! The stage index.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution values at x.
    real(real64), intent(inout), dimension(:) :: yn
        !!
    real(real64), intent(inout), dimension(:,:) :: fn
        !! The N-by-NSTAGES matrix containing the derivative values at each
        !! stage.
    integer(int32), intent(out) :: niter
        !! The number of Newton iterations performed.
    logical, intent(out) :: accept
        !! Returns true if the Newton iteration converged; else, false.
    real(real64), intent(out) :: resid
        !! The Euclidean norm of the residual.

    ! Local Variables
    integer(int32) :: j, n, maxiter
    real(real64) :: z, tol, alpha

    ! Initialization
    n = size(y)
    niter = 0
    accept = .false.
    call this%initialize_solver(n)
    maxiter = this%get_newton_iteration_limit()
    tol = this%get_newton_tolerance()
    alpha = this%get_model_coefficient(i, i)
    z = x + this%get_weight(i) * h

    ! Pre-Iteration Set-Up
    this%xi = 0.0d0
    if (i > 1) then
        do j = 1, i - 1
            this%xi = this%xi + this%get_model_coefficient(i,j) * fn(:,j)
        end do
        this%xi = y + h * this%xi
    end if
    yn = y
    call sys%fcn(z, yn, fn(:,i))

    ! Iteration
    do niter = 1, maxiter
        ! Compute the right-hand-side
        this%dy = this%xi + h * alpha * fn(:,i) - yn

        ! Solve the system
        call solve_lu(a, pvt, this%dy)

        ! Update the solution
        yn = yn + this%dy
        call sys%fcn(z, yn, fn(:,i))

        ! Check for convergence
        resid = norm2(this%dy)
        if (resid <= tol) then
            accept = .true.
            exit
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine dii_ode_solver(this, sys, x, iv, err)
    !! Solves the supplied system of ODE's.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array of independent variable values at which to return the 
        !! the solution to the ODE's.
    real(real64), intent(in), dimension(:) :: iv
        !! An array containing the initial values for each ODE.
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
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is 
        !!      defined.
        !!
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if there are less than 
        !!      2 values given in the independent variable array x.

    ! Local Variables
    logical :: dense, success, acceptNewton
    integer(int32) :: ii, i, j, n, neqn, flag, nsteps, nstages, niter, itrack, &
        maxiter
    integer(int32), allocatable, dimension(:) :: pvt
    real(real64) :: alpha, h, xo, xn, xmax, ei, eold, resid, rtrack
    real(real64), allocatable, dimension(:) :: f, y, yn, fn, yerr, yi, w
    real(real64), allocatable, dimension(:,:) :: k, jac, a
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    success = .true.

    ! Input Checking
    if (n < 2) then
        call report_array_size_error(errmgr, "dii_ode_solver", "x", 2, n)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "dii_ode_solver")
        return
    end if

    ! Additional Initialization
    neqn = size(iv)
    xo = x(1)
    xmax = x(n)
    j = 2
    eold = 1.0d-4
    dense = (n > 2)
    nsteps = this%get_step_limit()
    nstages = this%get_stage_count()
    alpha = this%get_model_coefficient(nstages, nstages)
    maxiter = this%get_newton_iteration_limit()

    ! Memory Allocations
    allocate( &
        f(neqn), &
        y(neqn), &
        yn(neqn), &
        fn(neqn), &
        yerr(neqn), &
        k(neqn, nstages),  &
        pvt(neqn), &
        jac(neqn, neqn), &
        a(neqn, neqn), &
        w(neqn), &
        stat = flag &
    )
    if (flag == 0 .and. dense) allocate(yi(neqn), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "dii_ode_solver", flag)
        return
    end if

    ! Estimate an initial step size
    !
    ! Outputs:
    ! - f: Value of the derivatives at xo
    ! - h: Initial step size estimate
    call this%estimate_inital_step_size(sys, xo, xmax, iv, f, h)
    k(:,1) = f
    
    ! Store the initial conditions
    call this%append_to_buffer(x(1), iv, errmgr)
    if (errmgr%has_error_occurred()) return
    y = iv

    ! Cycle until integration is complete
    outer : do i = 1, nsteps
        ! Perform any pre-step actions
        call this%pre_step_action(success, sys, h, xo, y, f, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Update the Jacobian, if necessary
        if (.not.this%get_is_matrix_current()) then
            call this%form_matrix(sys, alpha, h, xo, y, jac, a, pvt, errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Cycle over each stage
        itrack = 0
        rtrack = 0.0d0
        newton : do ii = 1, nstages
            ! Perform the Newton iteration
            call this%newton_iteration(sys, a, pvt, ii, h, xo, y, yn, k, &
                niter, acceptNewton, resid)
            itrack = max(niter, itrack)
            rtrack = max(resid, rtrack)
            if (.not.acceptNewton) exit newton
        end do newton
        
        ! If the Newton process failed, update the Jacobian and try again
        ! Also consider the case where convergence is difficult
        if (.not.acceptNewton) then
            call this%set_is_matrix_current(.false.)
            h = 0.5d0 * h   ! Half the step size as well
            if (abs(h) < abs(this%get_minimum_step_size())) then
                call report_step_size_too_small(errmgr, "dii_ode_solver", xo, h)
                return
            end if
            cycle outer
        end if

        ! Update the solution
        xn = xo + h
        do ii = 1, nstages
            if (ii == 1) then
                w = this%get_weight(ii) * k(:,ii)
            else
                w = w + this%get_weight(ii) * k(:,ii)
            end if
        end do
        yn = y + h * w

        ! Evaluate the error
        call this%compute_errors(sys, xo, y, f, xn, yn, fn, k, rtrack, itrack, &
            yerr)

        ! Compute a normalized error value.  A value of < 1 indicates success
        ei = this%compute_error_norm(y, yn, yerr)

        ! Determine the next step size
        h = this%estimate_next_step_size(ei, eold, h, xo, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Reject the step?
        success = ei <= 1.0d0
        if (.not.success) cycle outer ! We failed, try again with a new step size

        ! If we're here, the step has been successful.  Take any post-step
        ! action such as setting up interpolation routines, etc.
        call this%post_step_action(sys, dense, xo, xn, y, yn, f, fn, k)

        ! Do we need to interpolate for dense output, or can we just store
        ! values and move on
        if (dense) then
            ! Perform the interpolation as needed until a new step is required
            interp : do while (abs(x(j)) <= abs(xn))
                call this%interpolate(x(j), xo, y, f, xn, yn, fn, yi)
                call this%append_to_buffer(x(j), yi, errmgr)
                if (errmgr%has_error_occurred()) return
                j = j + 1
                if (j > n) exit interp
            end do interp
        else
            ! Store the values and move on
            call this%append_to_buffer(xn, yn, errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Are we done?
        if (abs(xn) >= abs(xmax)) then
            ! Deal with the case where output is only returned at the
            ! integration points and the solver oversteps the endpoint.
            if (abs(xn) > abs(xmax)) then
                ! Interpolate to get the solution at xmax
                call this%post_step_action(sys, .true., xo, xn, y, yn, f, fn, k)
                call this%interpolate(xmax, xo, y, f, xn, yn, fn, yi)
                call this%append_to_buffer(xmax, yi, errmgr)
                if (errmgr%has_error_occurred()) return
            end if
            
            ! We're done
            go to 100
        end if

        ! Do we need to limit the step size to not overshoot the terminal value?
        if (this%get_allow_overshoot() .and. abs(xn + h) > abs(xmax)) then
            ! Limit the step size
            h = xmax - xn
        end if

        ! Update parameters
        xo = xn
        y = yn

        ! Is this an FSAL integrator?  If so, we already have the derivative
        ! values.  If not, we need to recompute the derivative values for the
        ! next iteration
        if (this%get_is_fsal()) then
            ! We have the derivative estimate already
            f = fn
        else
            ! Update the derivative estimates
            call sys%fcn(xn, yn, f) ! write to f, not fn
        end if
        k(:,1) = f
    end do outer

    ! If we're here, the solver has run out of allowable steps
    call report_excessive_integration_steps(errmgr, "dii_ode_solver", nsteps, &
        xn)
    return

    ! End
100 continue
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine dii_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k)
    !! This is a pass-thru routine that has no function for DIRK integrators.
    class(diagonally_implicit_integrator), intent(inout) :: this
        !! The diagonally_implicit_integrator object.
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

    ! This is just a placeholder for this type of integrator as it is unused
    return
end subroutine

! ******************************************************************************
! 4TH ORDER IMPLICIT RUNGE-KUTTA
! ------------------------------------------------------------------------------
subroutine irk4_pre_step(this, prevs, sys, h, x, y, f, err)
    use diffeq_sdirk4_constants
    !! Placeholder routine for any pre-step actions.
    class(implicit_runge_kutta_4), intent(inout) :: this
        !! The implicit_runge_kutta_4 object.
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

    if (this%m_populated) return

    ! Populate the Butcher tableau
    this%m_a = 0.0d0

    this%m_a(2,1) = a21
    this%m_a(2,2) = a22

    this%m_a(3,1) = a31
    this%m_a(3,2) = a32
    this%m_a(3,3) = a33

    this%m_a(4,1) = a41
    this%m_a(4,2) = a42
    this%m_a(4,3) = a43
    this%m_a(4,4) = a44

    this%m_a(5,1) = a51
    this%m_a(5,2) = a52
    this%m_a(5,3) = a53
    this%m_a(5,4) = a54
    this%m_a(5,5) = a55

    this%m_a(6,1) = b1
    this%m_a(6,2) = b2
    this%m_a(6,3) = b3
    this%m_a(6,4) = b4
    this%m_a(6,5) = b5
    this%m_a(6,6) = b6

    ! B
    this%m_b(1) = b1
    this%m_b(2) = b2
    this%m_b(3) = b3
    this%m_b(4) = b4
    this%m_b(5) = b5
    this%m_b(6) = b6

    ! C
    this%m_c(1) = 0.0d0
    this%m_c(2) = c2
    this%m_c(3) = c3
    this%m_c(4) = c4
    this%m_c(5) = c5
    this%m_c(6) = c6

    ! Error Coefficients
    this%m_e(1) = b1a - b1
    this%m_e(2) = b2a - b2
    this%m_e(3) = b3a - b3
    this%m_e(4) = b4a - b4
    this%m_e(5) = b5a - b5
    this%m_e(6) = b6a - b6

    ! Interpolation Coefficients
    this%m_dc(1,1) = bs11
    this%m_dc(2,1) = bs21
    this%m_dc(3,1) = bs31
    this%m_dc(4,1) = bs41

    this%m_dc(1,2) = bs12
    this%m_dc(2,2) = bs22
    this%m_dc(3,2) = bs32
    this%m_dc(4,2) = bs42

    this%m_dc(1,3) = bs13
    this%m_dc(2,3) = bs23
    this%m_dc(3,3) = bs33
    this%m_dc(4,3) = bs43

    this%m_dc(1,4) = bs14
    this%m_dc(2,4) = bs24
    this%m_dc(3,4) = bs34
    this%m_dc(4,4) = bs44

    this%m_dc(1,5) = bs15
    this%m_dc(2,5) = bs25
    this%m_dc(3,5) = bs35
    this%m_dc(4,5) = bs45

    this%m_dc(1,6) = bs16
    this%m_dc(2,6) = bs26
    this%m_dc(3,6) = bs36
    this%m_dc(4,6) = bs46

    ! We're all set
    this%m_populated = .true.
end subroutine

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
    !! Gets the stage count for this integrator.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32) :: rst
        !! The stage count.
    rst = 6
end function

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
        !! The parameter.
    rst = this%m_a(i, j)
end function

! ------------------------------------------------------------------------------
pure function irk4_get_weight(this, i) result(rst)
    !! Gets the requested weighting factor from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    real(real64) :: rst
        !! The parameter.
    rst = this%m_b(i)
end function

! ------------------------------------------------------------------------------
pure function irk4_get_node(this, i) result(rst)
    !! Gets the requested node from the Butcher tableau.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    integer(int32), intent(in) :: i
        !! The row index.
    real(real64) :: rst
        !! The parameter.
    rst = this%m_c(i)
end function

! ------------------------------------------------------------------------------
subroutine irk4_eval_error(this, sys, x, y, f, xn, yn, fn, k, resid, &
    niter, yerr)
    !! Computes the error for each equation.
    class(implicit_runge_kutta_4), intent(in) :: this
        !! The implicit_runge_kutta_4 object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the derivatives at x.
    real(real64), intent(in) :: xn
        !! The next value of the independent variable.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution estimate at xn.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.
    real(real64), intent(in), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each
        !! stage.
    real(real64), intent(in) :: resid
        !! The norm of the residual from the Newton iteration process.
    integer(int32), intent(in) :: niter
        !! The number of iterations taken by the Newton process.
    real(real64), intent(out), dimension(:) :: yerr
        !! An N-element array where this routine will write the error
        !! estimates for each equation.

    ! Local Variables
    integer(int32) :: i, nstages
    real(real64) :: h

    ! Initialization
    nstages = this%get_stage_count()
    h = xn - x

    ! Process
    do i = 1, nstages
        if (i == 1) then
            yerr = this%m_e(i) * k(:,i)
        else
            yerr = yerr + this%m_e(i) * k(:,i)
        end if
    end do
    yerr = h * yerr
end subroutine

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

    ! Process
    if (allocated(this%m_k)) then
        if (size(this%m_k, 1) /= n .or. size(this%m_k, 2) /= nstages) then
            deallocate(this%m_k)
            allocate(this%m_k(n, nstages), source = k)
        else
            this%m_k = k
        end if
    else
        allocate(this%m_k(n, nstages), source = k)
    end if
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
    integer(int32) :: i,  j, norder, nstages, n
    real(real64) :: h, theta, bi
    real(real64), allocatable, dimension(:) :: yi

    ! Initialization
    n = size(yn)
    norder = this%get_order()
    nstages = this%get_stage_count()
    h = xn1 - xn
    theta = (x - xn) / h

    ! Process
    allocate(yi(n), source = 0.0d0)
    do i = 1, nstages
        bi = 0.0d0
        do j = 1, norder
            bi = bi + this%m_dc(j,i) * theta**j
        end do
        yi = yi + bi * this%m_k(:,i)
    end do
    y = yn + h * yi
end subroutine

! ! ------------------------------------------------------------------------------
! pure function irk4_estimate_error(this, y, yest, yerr) result(rst)
!     !! Computes the norm of the scaled error estimate.  A value less than one
!     !! indicates a successful step.  A value greater than one suggests that the
!     !! results do not meet the requested tolerances.
!     class(implicit_runge_kutta_4), intent(in) :: this
!         !! The implicit_runge_kutta_4 object.
!     real(real64), intent(in), dimension(:) :: y
!         !! The previously accepted solution array (N-element).
!     real(real64), intent(in), dimension(size(y)) :: yest
!         !! An N-element array containing the next solution point estimate.
!     real(real64), intent(in), dimension(size(y)) :: yerr
!         !! An N-element array containing the estimate of error for each
!         !! equation.
!     real(real64) :: rst
!         !! The norm of the scaled error.

!     ! Return 1 if the step size should be reduced, 0 if it's just right; else,
!     ! -1 if it should be increased
!     if (yerr(1) > 0.5d0 * this%get_newton_iteration_limit()) then
!         rst = 1.0d0
!     else if (yerr(1) < 0.15d0 * this%get_newton_iteration_limit()) then
!         rst = -1.0d0
!     else
!         rst = 0.0d0
!     end if
! end function

! ! ------------------------------------------------------------------------------
! function irk4_next_step(this, e, eold, h, x, err) result(rst)
!     !! Estimates the next step size based upon the current and previous error
!     !! estimates.
!     class(implicit_runge_kutta_4), intent(inout) :: this
!         !! The implicit_runge_kutta_4 object.
!     real(real64), intent(in) :: e
!         !! The norm of the current scaled error estimate.
!     real(real64), intent(inout) :: eold
!         !! The norm of the previous step's scaled error estimate.  On output,
!         !! this variable is updated.
!     real(real64), intent(in) :: h
!         !! The current step size.
!     real(real64), intent(in) :: x
!         !! The current independent variable value.
!     class(errors), intent(inout), optional, target :: err
!         !! An optional errors-based object that if provided 
!         !! can be used to retrieve information relating to any errors 
!         !! encountered during execution. If not provided, a default 
!         !! implementation of the errors class is used internally to 
!         !! provide error handling.  Possible errors and warning messages
!         !! that may be encountered are as follows.
!         !!
!         !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size
!         !!      becomes too small in magnitude.
!     real(real64) :: rst
!         !! The new step size estimate.

!     ! Local Variables
!     class(errors), pointer :: errmgr
!     type(errors), target :: deferr
!     real(real64) :: tol
    
!     ! Initialization
!     if (present(err)) then
!         errmgr => err
!     else
!         errmgr => deferr
!     end if
!     tol = sqrt(epsilon(tol))

!     ! This integrator bases step size adjustment on the efficiency of the
!     ! Newton solver.
!     if (abs(e - 1.0d0) <= tol) then
!         ! The solver had to work a bit.  It's worth reducing the step size
!         rst = 0.5d0 * h

!         ! Ensure the step size isn't too small
!         if (abs(rst) < abs(this%get_minimum_step_size())) then
!             call report_step_size_too_small(errmgr, "irk4_next_step", x, rst)
!             return
!         end if

!         ! Suggest a recomputation of the Jacobian
!         call this%set_is_matrix_current(.false.)
!     else if (abs(e + 1.0d0) <= tol) then
!         ! The solver had an easy go of things.  It's worth trying to increase
!         ! the step size
!         rst = 2.0d0 * h

!         ! Ensure the step size isn't too large
!         if (abs(rst) >= abs(this%get_maximum_step_size())) then
!             rst = sign(this%get_maximum_step_size(), h)
!         end if
!     else
!         ! Do nothing with the step size
!         rst = h
!     end if
! end function

! ------------------------------------------------------------------------------
end module