module diffeq_implicit_runge_kutta
    use iso_fortran_env
    use diffeq_runge_kutta
    use diffeq_errors
    use diffeq_base
    use linalg
    implicit none
    private
    public :: build_factored_newton_matrix_routine
    public :: implicit_rk_variable_integrator
    public :: sdirk_integrator
    public :: sdirk4_integrator

    type, abstract, extends(rk_variable_integrator) :: &
        implicit_rk_variable_integrator
        !! Defines an implicit Runge-Kutta variable-step integrator.
        real(real64), private :: m_hold = 0.0d0
            ! The most recent successful step size
        logical, private :: m_isJacCurrent = .false.
            ! Is the Jacobian current?
        logical, private :: m_usePI = .true.
            ! Use a PI step size (true) or a Gustafsson controller (false)
    contains
        procedure, public :: get_previous_step_size => irk_get_old_step
            !! Gets the most recent successful step size.
        procedure, public :: set_previous_step_size => irk_set_old_step
            !! Sets the most recent successful step size.
        procedure, public :: compute_next_step_size => &
            irk_compute_next_step_size
            !! Computes the next step size.
        procedure, public :: get_is_jacobian_current => irk_get_is_jac_current
            !! Gets a value determining if the Jacobian matrix estimate is
            !! current such that it does not need to be recomputed at this time.
        procedure, public :: set_is_jacobian_current => irk_set_is_jac_current
            !! Sets a value determining if the Jacobian matrix estimate is
            !! current such that it does not need to be recomputed at this time.
        procedure(build_factored_newton_matrix_routine), public, deferred :: &
            build_factored_newton_matrix
            !! Builds the matrix of the form \( X = f I - J \) or 
            !! \( X = f M - J \) if a mass matrix is defined, and then computes
            !! its LU factorization.  The Jacobian and mass matrices are 
            !! evaluated as part of this process, if necessary.
        procedure, public :: get_use_pi_controller => irk_get_use_pi_controller
            !! Gets a parameter determining if a PI step size controller
            !! or a Gustafsson step size controller should be used. The default
            !! is to use a PI step size controller.
        procedure, public :: set_use_pi_controller => irk_set_use_pi_controller
            !! Sets a parameter determining if a PI step size controller
            !! or a Gustafsson step size controller should be used. The default
            !! is to use a PI step size controller.
        procedure, public :: estimate_first_step_size => irk_estimate_first_step
            !! Computes an estimate to the first step size based upon the
            !! initial function values.
    end type

    interface
        subroutine build_factored_newton_matrix_routine(this, sys, h, x, y, err)
            use iso_fortran_env, only : real64
            use ferror, only : errors
            import implicit_rk_variable_integrator
            import ode_container
            class(implicit_rk_variable_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, abstract, extends(implicit_rk_variable_integrator) :: sdirk_integrator
        !! Defines a base structure for singly diagonally implicit 
        !! Runge-Kutta integrators.
        real(real64), private, allocatable, dimension(:,:) :: m_jac
            ! Jacobian matrix workspace
        real(real64), private, allocatable, dimension(:,:) :: m_mass
            ! Mass matrix workspace
        real(real64), private, allocatable, dimension(:,:) :: m_mtx
            ! System matrix workspace
        integer(int32), private, allocatable, dimension(:) :: m_pvt
            ! LU pivot tracking workspace
        real(real64), private, allocatable, dimension(:) :: m_w
            ! NEQN-element workspace array
        real(real64), private, allocatable, dimension(:) :: m_dy
            ! NEQN-element workspace array
        integer(int32), private :: m_maxNewtonIter = 7
            ! Allowable number of Newton iterations
        real(real64), private :: m_newtontol = 1.0d-6
            ! Newton iteration tolerance
    contains
        procedure, public :: initialize => sdirk_alloc_workspace
            !! Initializes the integrator.
        procedure, public :: build_newton_matrix => sdirk_build_matrix
            !! Builds the system matrix of the form \( X = f I - J \)
            !! or \( X = f M - J \) if a mass matrix is defined.
        procedure, public :: build_factored_newton_matrix => &
            sdirk_build_factored_matrix
            !! Builds the matrix of the form \( X = f I - J \) or 
            !! \( X = f M - J \) if a mass matrix is defined, and then computes
            !! its LU factorization.  The Jacobian and mass matrices are 
            !! evaluated as part of this process, if necessary.
        procedure, public :: get_max_newton_iteration_count => &
            sdirk_get_max_newton_iter
            !! Gets the maximum allowed number of Newton iterations.
        procedure, public :: set_max_newton_iteration_count => &
            sdirk_set_max_newton_iter
            !! Sets the maximum allowed number of Newton iterations.
        procedure, public :: get_newton_tolerance => sdirk_get_newton_tol
            !! Gets the tolerance used to check for convergence of the
            !! Newton iterations.
        procedure, public :: set_newton_tolerance => sdirk_set_newton_tol
            !! Sets the tolerance used to check for convergence of the
            !! Newton iterations.
        procedure, public :: attempt_step => sdirk_attempt_step
            !! Attempts a single integration step.
        procedure, public :: solve_newton_stage => sdirk_solve_newton
            !! Solves the Newton iteration problem for the i-th stage.
    end type

! ------------------------------------------------------------------------------
    type, extends(sdirk_integrator) :: sdirk4_integrator
        !! Defines a singly diagonally implicit 4th order Runge-Kutta 
        !! integrator suitable for integrating stiff systems of differential 
        !! equations.
        logical, private :: m_modelDefined = .false.
        real(real64), private, allocatable, dimension(:,:) :: m_dc
            ! Interpolation coefficients for dense output
    contains
        procedure, public :: get_order => sd4_get_order
            !! Returns the order of the integrator.
        procedure, public :: get_stage_count => sd4_get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure, public :: define_model => sd4_define_model
            !! Defines (initializes) the model parameters.
        procedure, public :: is_fsal => sd4_is_fsal
            !! Determines if the integrator is an FSAL (first same as last)
            !! integrator.
        procedure, public :: interpolate => sd4_interp
            !! Provides interpolation between integration points allowing for
            !! dense output.
        procedure, public :: set_up_interpolation => sd4_set_up_interp
            !! Sets up the interpolation polynomial.
        procedure, public :: initialize => sd4_alloc_workspace
            !! Initializes the integrator.
    end type

contains
! ------------------------------------------------------------------------------
pure function irk_get_old_step(this) result(rst)
    !! Gets the most recent successful step size.
    class(implicit_rk_variable_integrator), intent(in) :: this
        !! The implicit_rk_variable_integrator object.
    real(real64) :: rst
        !! The step size.
    rst = this%m_hold
end function

! ---------------------
subroutine irk_set_old_step(this, x)
    !! Sets the most recent successful step size.
    class(implicit_rk_variable_integrator), intent(inout) :: this
        !! The implicit_rk_variable_integrator object.
    real(real64), intent(in) :: x
        !! The step size.
    this%m_hold = x
end subroutine

! ------------------------------------------------------------------------------
function irk_compute_next_step_size(this, hn, en, enm1) result(rst)
    !! Computes the next step size.
    class(implicit_rk_variable_integrator), intent(inout) :: this
        !! The implicit_rk_variable_integrator object.
    real(real64), intent(in) :: hn
        !! The current step size.
    real(real64), intent(in) :: en
        !! The norm of the error for the current step size.
    real(real64), intent(in) :: enm1
        !! The norm of the error from the previous step size.
    real(real64) :: rst
        !! The new step size.

    ! Parameters
    real(real64), parameter :: fac1 = 5.0d0
    real(real64), parameter :: fac2 = 1.0d0 / 8.0d0

    ! Local Variables
    real(real64) :: fac, safe, hnew, facpred, e, hold, hmax

    ! Process
    if (this%get_use_pi_controller()) then
        rst = pi_controller(this%get_alpha(), this%get_beta(), &
            this%get_order(), hn, en, enm1, this%get_safety_factor(), &
            this%get_max_step_size())
    else
        ! Initialization
        safe = this%get_safety_factor()
        e = 1.0d0 / real(this%get_order(), real64)
        fac = max(fac2, min(fac1, (en**2 / enm1)**e / safe))
        rst = hn / fac
        hold = this%get_previous_step_size()
        hmax = this%get_max_step_size()

        ! Process
        if (en <= 1.0d0) then
            if (.not.this%get_is_first_step()) then
                facpred = (hold / hn) * (en**2 / enm1)**e / safe
                facpred = max(fac2, min(fac1, facpred))
                fac = max(fac, facpred)
                rst = hn / fac
            end if
            call this%set_is_first_step(.false.)
            call this%set_previous_step_size(hn)
            call this%set_previous_error_norm(max(1.0d-2, en))
        end if
        if (abs(rst) > abs(hmax)) rst = sign(hmax, rst)
    end if
end function

! ------------------------------------------------------------------------------
pure function irk_get_is_jac_current(this) result(rst)
    !! Gets a value determining if the Jacobian matrix estimate is current such 
    !! that it does not need to be recomputed at this time.
    class(implicit_rk_variable_integrator), intent(in) :: this
        !! The implicit_rk_variable_integrator object.
    logical :: rst
        !! True if the Jacobian matrix is current; else, false.
    rst = this%m_isJacCurrent
end function

! ---------------------
subroutine irk_set_is_jac_current(this, x)
    !! Sets a value determining if the Jacobian matrix estimate is current such 
    !! that it does not need to be recomputed at this time.
    class(implicit_rk_variable_integrator), intent(inout) :: this
        !! The implicit_rk_variable_integrator object.
    logical, intent(in) :: x
        !! True if the Jacobian matrix is current; else, false.
    this%m_isJacCurrent = x
end subroutine

! ------------------------------------------------------------------------------
pure function irk_get_use_pi_controller(this) result(rst)
    !! Gets a parameter determining if a PI step size controller or a Gustafsson
    !! step size controller should be used. The default is to use a PI step 
    !! size controller.
    class(implicit_rk_variable_integrator), intent(in) :: this
        !! The implicit_rk_variable_integrator object.
    logical :: rst
        !! True to use a PI controller; else, false to use a Gustafsson 
        !! controller.
    rst = this%m_usePI
end function

! ---------------------
subroutine irk_set_use_pi_controller(this, x)
    !! Sets a parameter determining if a PI step size controller or a Gustafsson
    !! step size controller should be used. The default is to use a PI step 
    !! size controller.
    class(implicit_rk_variable_integrator), intent(inout) :: this
        !! The implicit_rk_variable_integrator object.
    logical, intent(in) :: x
        !! True to use a PI controller; else, false to use a Gustafsson 
        !! controller.
    this%m_usePI = x
end subroutine

! ------------------------------------------------------------------------------
pure function irk_estimate_first_step(this, xo, xf, yo, fo) result(rst)
    !! Computes an estimate to the first step size based upon the initial 
    !! function values.
    class(implicit_rk_variable_integrator), intent(in) :: this
        !! The implicit_rk_variable_integrator object.
    real(real64), intent(in) :: xo
        !! The initial value of the independent variable.
    real(real64), intent(in) :: xf
        !! The final value of the independent variable.
    real(real64), intent(in), dimension(:) :: yo
        !! An N-element array containing the initial values.
    real(real64), intent(in), dimension(:) :: fo
        !! An N-element array containing the initial function values.
    real(real64) :: rst
        !! An estimate on the initial step size.

    ! Local Variables
    real(real64) :: h1, h2

    ! Process
    h1 = 1.0d-5 * (xf - xo)
    h2 = this%get_max_step_size()
    rst = sign(min(abs(h1), abs(h2)), h1)
end function

! ******************************************************************************
! SDIRK INTEGRATOR - BASE TYPE
! ------------------------------------------------------------------------------
subroutine sdirk_alloc_workspace(this, neqn, err)
    !! Initializes the integrator.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.  Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.

    ! Local Variables
    integer(int32) :: flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Allocations
    if (allocated(this%m_jac)) then
        if (size(this%m_jac, 1) /= neqn .or. size(this%m_jac, 2) /= neqn) then
            deallocate(this%m_jac)
            allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_mass)) then
        if (size(this%m_mass, 1) /= neqn .or. size(this%m_mass, 2) /= neqn) then
            deallocate(this%m_mass)
            allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_mtx)) then
        if (size(this%m_mtx, 1) /= neqn .or. size(this%m_mtx, 2) /= neqn) then
            deallocate(this%m_mtx)
            allocate(this%m_mtx(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_mtx(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_pvt)) then
        if (size(this%m_pvt) /= neqn) then
            deallocate(this%m_pvt)
            allocate(this%m_pvt(neqn), stat = flag, source = 0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_pvt(neqn), stat = flag, source = 0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_w)) then
        if (size(this%m_w) /= neqn) then
            deallocate(this%m_w)
            allocate(this%m_w(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_w(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_dy)) then
        if (size(this%m_dy) /= neqn) then
            deallocate(this%m_dy)
            allocate(this%m_dy(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_dy(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! Ensure derivative storage is set up
    call this%initialize_derivative_storage(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Set up the model
    call this%initialize_model_storage(errmgr)
    if (errmgr%has_error_occurred()) return

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "sdirk_alloc_workspace", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine sdirk_build_factored_matrix(this, sys, h, x, y, err)
    !! Builds the matrix of the form \( X = f I - J \) or 
    !! \( X = f M - J \) if a mass matrix is defined, and then computes
    !! its LU factorization.  The Jacobian and mass matrices are 
    !! evaluated as part of this process, if necessary.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An array containing the current values of the dependent variables.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.  Possible errors and warning messages that may be 
        !! encountered are as follows.

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    logical :: change
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    change = .false.

    ! Compute the Jacobian, if an update is needed.
    if (.not.this%get_is_jacobian_current() .or. this%get_is_first_step()) then
        ! Recompute the Jacobian
        call sys%compute_jacobian(x, y, this%m_jac, errmgr)
        if (errmgr%has_error_occurred()) return
        call this%set_is_jacobian_current(.true.)
        change = .true.
    end if

    ! Compute the mass matrix
    if (associated(sys%mass_matrix) .and. &
        (sys%get_is_mass_matrix_dependent() .or. this%get_is_first_step())) &
    then
        call sys%mass_matrix(x, y, this%m_mass)
        change = .true.
    end if

    ! Compute the system matrix
    if (associated(sys%mass_matrix)) then
        call this%build_newton_matrix(h, this%m_jac, this%m_mtx, this%m_mass)
    else
        call this%build_newton_matrix(h, this%m_jac, this%m_mtx)
    end if

    ! Compute the LU factorization of the system matrix if anything has changed
    if (change) then
        call lu_factor(this%m_mtx, this%m_pvt, errmgr)
        if (errmgr%has_error_occurred()) return
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine sdirk_build_matrix(this, h, jac, x, m, err)
    !! Builds the system matrix of the form \( X = f I - J \)
    !! or \( X = f M - J \) if a mass matrix is defined.
    class(sdirk_integrator), intent(in) :: this
        !! The sdirk_integrator object.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in), dimension(:,:) :: jac
        !! The current NEQN-by-NEQN Jacobian matrix.
    real(real64), intent(out), dimension(:,:) :: x
        !! An NEQN-by-NEQN matrix where the output will be written.
    real(real64), intent(in), dimension(:,:), optional :: m
        !! An optional NEQN-by-NEQN mass matrix.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.  Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !!  - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if any of the matrices are not
        !!      sized correctly.

    ! Local Variables
    integer(int32) :: i, neqn, ns
    real(real64) :: fac
    class(errors), pointer :: errmgr
    type(errors), target :: deferr

    ! Initialization
    neqn = size(this%m_jac, 1)
    ns = this%get_stage_count()
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Checking
    if (size(jac, 1) /= neqn .or. size(jac, 2) /= neqn) then
        call report_matrix_size_error(errmgr, "sdirk_build_matrix", "jac", &
            neqn, neqn, size(jac, 1), size(jac, 2))
        return
    end if
    if (size(x, 1) /= neqn .or. size(x, 2) /= neqn) then
        call report_matrix_size_error(errmgr, "sdirk_build_matrix", "x", &
            neqn, neqn, size(x, 1), size(x, 2))
        return
    end if
    if (present(m)) then
        if (size(m, 1) /= neqn .or. size(m, 2) /= neqn) then
            call report_matrix_size_error(errmgr, "sdirk_build_matrix", "m", &
                neqn, neqn, size(m, 1), size(m, 2))
            return
        end if
    end if

    ! Process
    fac = this%a(ns, ns) * h
    if (present(m)) then
        x = m - fac * jac
    else
        x = -fac * jac
        do i = 1, neqn
            x(i,i) = x(i,i) + 1.0d0
        end do
    end if

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine sdirk_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
    fprev, err)
    !! Attempts a single integration step.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the dependent 
        !! variables.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where the values of the dependent variables at 
        !! x + h will be written.
    real(real64), intent(out), dimension(:) :: en
    real(real64), intent(in), optional, dimension(:) :: xprev
        !! An optional M-element array containing the previous M values
        !! of the independent variable where M is the order of the 
        !! method.  This is typically only used for multi-step methods.
        !! In single-step methods, this parameter is typically not
        !! needed.
    real(real64), intent(in), optional, dimension(:,:) :: yprev
        !! An optional M-by-N matrix containing the previous M arrays of
        !! dependent variables, where M is the order of the method.  As
        !! with xprev, this parameter is typically used for multi-step
        !! methods.  In single-step methods, this parameter is typically
        !! not needed.
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
        !! An optional M-by-N matrix containing the previous M arrays of
        !! ODE (function) values.  As with xprev and yprev, M is the
        !! order of the method, and this parameter is typically used for
        !! multi-step methods.  In single-step methods, this parameter 
        !! is typically not needed.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.

    ! Local Variables
    logical :: accept
    integer(int32) :: i, j, neqn, nstages, niter, maxiter, itertracking
    real(real64) :: z, tol, disp, val, eval
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    nstages = this%get_stage_count()
    maxiter = this%get_max_newton_iteration_count()
    accept = .false.
    tol = this%get_newton_tolerance()
    itertracking = 0

    ! Initialize, if necessary
    if (.not.allocated(this%m_jac)) then
        call this%initialize(neqn, errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! Ensure the Jacobian is up to date prior to the step
    if (.not.this%get_is_jacobian_current()) then
        call this%build_factored_newton_matrix(sys, h, x, y, errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! Process
    if (.not.this%is_fsal() .or. this%get_is_first_step()) then
        ! On FSAL integrators, we only need to make this call on the first step
        ! as the integrator uses the last evaluation from the previous step
        ! as this step.  On non-FSAL integrators we always need to compute an
        ! updated first step.
        call sys%ode(x, y, this%f(:,1))
    end if

    ! Cycle over each stage and solve the Newton problem
    do i = 1, nstages
        ! Attempt to solve the Newton problem
        call this%solve_newton_stage(sys, i, h, x, y, yn, accept, niter)
        itertracking = max(niter, itertracking)
        if (.not.accept) exit
    end do

    ! Do we need to update the Jacobian?
    if (.not.accept) then
        ! We couldn't converge - force a Jacobian update
        call this%set_is_jacobian_current(.false.)

        ! TO DO: 
        ! Figure out a way to estimate a new step size 
    end if

    ! Update the solution estimate and error estimate
    do i = 1, nstages
        if (i == 1) then
            yn = this%b(i) * this%f(:,i)
            en = this%e(i) * this%f(:,i)
        else
            yn = yn + this%b(i) * this%f(:,i)
            en = en + this%e(i) * this%f(:,i)
        end if
    end do
    yn = y + h * yn
    en = h * en

    ! Only update the Jacobian if the Newton iteration encountered difficulty
    if (itertracking > maxiter / 2) then
        call this%set_is_jacobian_current(.false.)
    end if

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
pure function sdirk_get_max_newton_iter(this) result(rst)
    !! Gets the maximum allowed number of Newton iterations.
    class(sdirk_integrator), intent(in) :: this
        !! The sdirk_integrator object.
    integer(int32) :: rst
        !! The iteration limit.
    rst = this%m_maxNewtonIter
end function

! --------------------
subroutine sdirk_set_max_newton_iter(this, x)
    !! Sets the maximum allowed number of Newton iterations.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    integer(int32), intent(in) :: x
        !! The iteration limit.
    this%m_maxNewtonIter = x
end subroutine

! ------------------------------------------------------------------------------
pure function sdirk_get_newton_tol(this) result(rst)
    !! Gets the tolerance used to check for convergence of the Newton 
    !! iterations.
    class(sdirk_integrator), intent(in) :: this
        !! The sdirk_integrator object.
    real(real64) :: rst
        !! The tolerance.
    rst = this%m_newtontol
end function

! --------------------
subroutine sdirk_set_newton_tol(this, x)
    !! Sets the tolerance used to check for convergence of the Newton 
    !! iterations.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    real(real64), intent(in) :: x
        !! The tolerance.
    this%m_newtontol = x
end subroutine

! ------------------------------------------------------------------------------
subroutine sdirk_solve_newton(this, sys, i, h, x, y, yw, accept, niter)
    !! Solves the Newton iteration problem for the i-th stage.
    class(sdirk_integrator), intent(inout) :: this
        !! The sdirk_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    integer(int32), intent(in) :: i
        !! The current stage number.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the dependent
        !! variables.
    real(real64), intent(out), dimension(:) :: yw
        !! An N-element workspace array.
    logical, intent(out) :: accept
        !! Returns true if the Newton iteration reached convergence; else, 
        !! false if the iteration did not converge.
    integer(int32), intent(out) :: niter
        !! The number of iterations performed.

    ! Local Variables
    integer(int32) :: j, maxiter
    real(real64) :: z, disp, tol, alpha

    ! Initialization
    accept = .false.
    z = x + this%c(i) * h
    tol = this%get_newton_tolerance()
    maxiter = this%get_max_newton_iteration_count()
    yw = y

    ! Process
    this%m_w = 0.0d0
    do j = 1, i - 1
        alpha = this%a(i,j)
        this%m_w = this%m_w + alpha * this%f(:,i)
    end do
    this%m_w = y + h * this%m_w
    call sys%ode(z, this%m_w, this%f(:,i))

    do niter = 1, maxiter
        ! Compute the right-hand-side
        this%m_dy = this%m_w + h * this%a(i,i) * this%f(:,i) - yw

        ! Solve the system
        call solve_lu(this%m_mtx, this%m_pvt, this%m_dy)

        ! Update the solution
        yw = yw + this%m_dy

        ! Update the function evaluation
        call sys%ode(z, yw, this%f(:,i))

        ! Check for convergence
        disp = norm2(this%m_dy)
        if (disp < tol) then
            accept = .true.
            exit
        end if
    end do
end subroutine

! ******************************************************************************
! SINGLY DIAGONALLY IMPLICIT 4TH ORDER RUNGE-KUTTA INTEGRATOR
! ------------------------------------------------------------------------------
subroutine sd4_alloc_workspace(this, neqn, err)
    use diffeq_sdirk4_constants
    !! Initializes the integrator.
    class(sdirk4_integrator), intent(inout) :: this
        !! The sdirk4_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.  Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !! - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.

    ! Local Variables
    integer(int32) :: norder, nstages, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    norder = this%get_order()
    nstages = this%get_stage_count()

    ! Process
    if (allocated(this%m_dc)) then
        if (size(this%m_dc, 1) /= norder .or. size(this%m_dc, 2) /= nstages) &
        then
            deallocate(this%m_dc)
            allocate(this%m_dc(norder, nstages), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_dc(norder, nstages), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! Populate the coefficient matrix
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

    ! Call the base routine
    call sdirk_alloc_workspace(this, neqn, errmgr)

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "sd4_alloc_workspace", flag)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
pure function sd4_get_order(this) result(rst)
    !! Returns the order of the integrator.
    class(sdirk4_integrator), intent(in) :: this
        !! The sdirk4_integrator object.
    integer(int32) :: rst
        !! The order of the integrator, 4 in this case.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function sd4_get_stage_count(this) result(rst)
    !! Gets the number of stages used by the integrator.
    class(sdirk4_integrator), intent(in) :: this
        !! The sdirk4_integrator object.
    integer(int32) :: rst
        !! The number of stages, 6 in this case.
    rst = 6
end function

! ------------------------------------------------------------------------------
subroutine sd4_define_model(this)
    use diffeq_sdirk4_constants
    !! Defines (initializes) the model parameters.
    class(sdirk4_integrator), intent(inout) :: this
        !! The sdirk4_integrator object.

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%a = 0.0d0

    this%a(2,1) = a21
    this%a(2,2) = a22

    this%a(3,1) = a31
    this%a(3,2) = a32
    this%a(3,3) = a33

    this%a(4,1) = a41
    this%a(4,2) = a42
    this%a(4,3) = a43
    this%a(4,4) = a44

    this%a(5,1) = a51
    this%a(5,2) = a52
    this%a(5,3) = a53
    this%a(5,4) = a54
    this%a(5,5) = a55

    this%a(6,1) = b1
    this%a(6,2) = b2
    this%a(6,3) = b3
    this%a(6,4) = b4
    this%a(6,5) = b5
    this%a(6,6) = b6

    ! B
    this%b(1) = b1
    this%b(2) = b2
    this%b(3) = b3
    this%b(4) = b4
    this%b(5) = b5
    this%b(6) = b6

    ! C
    this%c(1) = 0.0d0
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4
    this%c(5) = c5
    this%c(6) = c6

    ! E
    this%e(1) = b1a - b1
    this%e(2) = b2a - b2
    this%e(3) = b3a - b3
    this%e(4) = b4a - b4
    this%e(5) = b5a - b5
    this%e(6) = b6a - b6

    ! Update definition status
    this%m_modelDefined = .true.
end subroutine

! ------------------------------------------------------------------------------
pure function sd4_is_fsal(this) result(rst)
    !! Determines if the integrator is an FSAL (first same as last) integrator.
    class(sdirk4_integrator), intent(in) :: this
        !! The sdirk4_integrator object.
    logical :: rst
        !! Returns true as this is an FSAL integrator.
    rst = .true.
end function

! ------------------------------------------------------------------------------
subroutine sd4_set_up_interp(this, x, xn, y, yn, k)
    !! Sets up the interpolation polynomial.
    class(sdirk4_integrator), intent(inout) :: this
        !! The sdirk4_integrator object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in) :: xn
        !! The value of the independent variable at the next step.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution values.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution values at the next step.
    real(real64), intent(in), dimension(:,:) :: k
        !! An N-by-M matrix containing the intermediate step function outputs 
        !! where M is the number of stages of the integrator.

    ! No set-up actions required
end subroutine

! ------------------------------------------------------------------------------
subroutine sd4_interp(this, xprev, yprev, xnew, x, y, err)
    !! Provides interpolation between integration points allowing for dense 
    !! output.
    class(sdirk4_integrator), intent(in) :: this
        !! The sdirk4_integrator object.
    real(real64), intent(in) :: xprev
        !! The previous value of the independent variable.
    real(real64), intent(in), dimension(:) :: yprev
        !! An N-element array containing the values of the dependent variables 
        !! at xprev.
    real(real64), intent(in) :: xnew
        !! The updated value of the independent variable.
    real(real64), intent(in) :: x
        !! The value at which to perform the interpolation.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array containing the interpolated values for each 
        !! equation.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling.

    ! Local Variables
    integer(int32) :: i, j, norder, nstages, neqn
    real(real64) :: h, theta, bi
    real(real64), allocatable, dimension(:) :: yn

    ! Initialization
    neqn = size(yprev)
    norder = this%get_order()
    nstages = this%get_stage_count()
    h = xnew - xprev
    theta = (x - xprev) / h

    ! Process
    allocate(yn(neqn), source = 0.0d0)
    do i = 1, nstages
        bi = 0.0d0
        do j = 1, norder
            bi = bi + this%m_dc(j,i) * theta**j
        end do

        yn = yn + bi * this%f(:,i)
    end do
    y = yprev + h * yn
end subroutine

! ------------------------------------------------------------------------------
end module