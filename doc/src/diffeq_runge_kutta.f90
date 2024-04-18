module diffeq_runge_kutta
    use iso_fortran_env
    use diffeq_variable_singlestep
    use diffeq_errors
    use diffeq_base
    implicit none
    private
    public :: pi_controller
    public :: rkv_get_matrix_parameter
    public :: rkv_get_array_parameter
    public :: rkv_get_boolean_parameter
    public :: rkv_get_integer_parameter
    public :: rkv_action
    public :: rkv_set_up_interp
    public :: rk_variable_integrator
    public :: dprk45_integrator
    public :: bsrk32_integrator

    !> @brief Defines a variable-step, Runge-Kutta integrator.
    type, abstract, extends(variable_singlestep_integrator) :: &
        rk_variable_integrator
        ! Workspace matrix (NEQN -by- STAGE COUNT)
        real(real64), public, allocatable, dimension(:,:) :: f
            !! An NEQN-by-NSTAGE matrix containing the function evaluations
            !! (derivatives) at each of the stages of evaluation.
        real(real64), private, allocatable, dimension(:) :: m_ywork
            ! Workspace array (NEQN)
        logical, private :: m_firstStep = .true.
            ! A flag determining if this is the first accepted step (use for FSAL)
        real(real64), private :: m_alpha = 0.7d0
            ! Step-size PI control parameters
        real(real64), private :: m_beta = 0.4d-1
            ! Step-size PI control parameters
        real(real64), public, allocatable, dimension(:,:) :: a
            !! The NSTAGE-by-NSTAGE method factor matrix A.
        real(real64), public, allocatable, dimension(:) :: b
            !! The NSTAGE-element quadrature weight array B.
        real(real64), public, allocatable, dimension(:) :: c
            !! The NSTAGE-element position factor array C.
        real(real64), public, allocatable, dimension(:) :: e
            !! The NSTAGE-element error estimate factor array E.
    contains
        procedure, public :: initialize => rkv_alloc_workspace
            !! Initializes the integrator.
        procedure(rkv_get_boolean_parameter), deferred, public :: is_fsal
            !! Gets a value determining if this is a FSAL (first, same as last)
            !! integrator.
        procedure(rkv_get_integer_parameter), deferred, public :: &
            get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure(rkv_action), deferred, public :: define_model
            !! Defines (initializes) the model parameters.
        procedure, public :: reset => rkv_reset
            !! Resets the integrator to its initial state.
        procedure, public :: attempt_step => rkv_attempt_step
            !! Attempts a single integration step.
        procedure, public :: on_successful_step => rkv_on_successful_step
            !! Perform necessary actions on completion of a successful step.
        procedure(rkv_set_up_interp), public, deferred :: set_up_interpolation
            !! Sets up the interpolation polynomial.
        procedure, public :: get_alpha => rkv_get_alpha
            !! Gets the \( \alpha \) control parameter in the PI controller
            !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
            !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
            !! factor, and \f$ k \) is the order of the integration method.
        procedure, public :: set_alpha => rkv_set_alpha
            !! Sets the \( \alpha \) control parameter in the PI controller
            !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
            !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
            !! factor, and \f$ k \) is the order of the integration method.
        procedure, public :: get_beta => rkv_get_beta
            !! Gets the \( \beta \) control parameter in the PI controller
            !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
            !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
            !! factor, and \f$ k \) is the order of the integration method.
        procedure, public :: set_beta => rkv_set_beta
            !! Sets the \( \beta \) control parameter in the PI controller
            !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
            !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
            !! factor, and \f$ k \) is the order of the integration method.
        procedure, public :: compute_next_step_size => rkv_next_step
            !! Computes the next step size.
        procedure, public :: get_is_first_step => rkv_get_is_first_step
            !! Gets a value determining if the integrator is set up to take its 
            !! first integration step.
        procedure, public :: set_is_first_step => rkv_set_is_first_step
            !! Sets a value determining if the integrator is set up to take its 
            !! first integration step.
    end type

    interface
        pure function rkv_get_matrix_parameter(this, i, j) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32), intent(in) :: i, j
            real(real64) :: rst
        end function

        pure function rkv_get_array_parameter(this, i) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32), intent(in) :: i
            real(real64) :: rst
        end function

        pure function rkv_get_boolean_parameter(this) result(rst)
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            logical :: rst
        end function

        pure function rkv_get_integer_parameter(this) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        subroutine rkv_action(this)
            !! Defines an action to undertake on a rk_variable_integrator 
            !! object.
            import rk_variable_integrator
            class(rk_variable_integrator), intent(inout) :: this
                !! The rk_variable_integrator object.
        end subroutine

        subroutine rkv_set_up_interp(this, x, xn, y, yn, k)
            !! Sets up interpolation for the rk_variable_integrator object.
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(inout) :: this
                !! The rk_variable_integrator object.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in) :: xn
                !! The value of the independent variable at x + h.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the values of the dependent
                !! variables at x.
            real(real64), intent(in), dimension(:) :: yn
                !! An N-element array containing the values of the dependent
                !! variables at x + h.
            real(real64), intent(in), dimension(:,:) :: k
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, extends(rk_variable_integrator) :: dprk45_integrator
        !! Defines the classical Dormand-Prince 4th/5th order integrator.
        logical, private :: m_modelDefined = .false.
        real(real64), private, allocatable, dimension(:,:) :: m_dprk45work
    contains
        procedure, public :: define_model => dprk45_define_model
            !! Defines (initializes) the model parameters.
        procedure, public :: is_fsal => dprk45_is_fsal
            !! Determines if the integrator is an FSAL (first same as last)
            !! integrator (e.g. the 4th/5th order Dormand-Prince integrator).
        procedure, public :: get_stage_count => dprk45_get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure, public :: get_order => dprk45_get_order
            !! Returns the order of the integrator.
        procedure, public :: interpolate => dprk45_interp
            !! Provides interpolation between integration points allowing
            !! for dense output.
        procedure, public :: set_up_interpolation => dprk45_set_up_interp
            !! Sets up the interpolation polynomial.
    end type

! ------------------------------------------------------------------------------
    type, extends(rk_variable_integrator) :: bsrk32_integrator
        !! Defines the Bogacki-Shampine 3rd/2nd order integrator.
        logical, private :: m_modelDefined = .false.
        real(real64), private, allocatable, dimension(:,:) :: m_bsrk23work
    contains
        procedure, public :: define_model => bsrk32_define_model
            !! Defines (initializes) the model parameters.
        procedure, public :: is_fsal => bsrk32_is_fsal
            !! Determines if the integrator is an FSAL (first same as last)
            !! integrator.
        procedure, public :: get_stage_count => bsrk32_get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure, public :: get_order => bsrk32_get_order
            !! Returns the order of the integrator.
        procedure, public :: interpolate => bsrk32_interp
            !! Provides interpolation between integration points allowing for
            !! dense output.
        procedure, public :: set_up_interpolation => bsrk32_set_up_interp
            !! Sets up the interpolation polynomial.
    end type

contains
! ------------------------------------------------------------------------------
pure function pi_controller(alpha, beta, order, hn, en, enm1, fs, maxstep) &
    result(rst)
    !! Computes the next step size using a PI type controller such that
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
    !! e_{n-1}^{\beta} \) where \( f \) is the safety factor and \( k \)
    !! is the order of the integration method.
    real(real64), intent(in) :: alpha
        !! The \( \alpha \) control parameter.
    real(real64), intent(in) :: beta
        !! The \( \beta \) control parameter.
    integer(int32), intent(in) :: order
        !! The order of the integrator.
    real(real64), intent(in) :: hn
        !! The current step size.
    real(real64), intent(in) :: en
        !! The norm of the error for the current step size.
    real(real64), intent(in) :: enm1
        !! The norm of the error from the previous step size.
    real(real64), intent(in) :: fs
        !! A safety factor to place on the growth of the step size.
    real(real64), intent(in) :: maxstep
        !! A cap on the size of the maximum step.
    real(real64) :: rst
        !! The new step size.

    ! Local Variables
    real(real64) :: hest

    ! Process
    hest = fs * hn * (1.0d0 / en)**(1.0d0 / order)
    rst = hest * (en**alpha) * (enm1**beta)

    if (abs(rst) > maxstep) rst = sign(maxstep, rst)
    if (rst / hn > 2.0d0) rst = 2.0d0 * hn
end function

! ******************************************************************************
! ABSTRACT VARIABLE-STEP RUNGE-KUTTA
! ------------------------------------------------------------------------------
subroutine rkv_alloc_workspace(this, neqn, err)
    !! Initializes the integrator.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations being integrated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.

    ! Local Variables
    integer(int32) :: m, n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Process
    m = neqn
    n = this%get_stage_count()
    if (allocated(this%f)) then
        if (size(this%f, 1) /= m .or. size(this%f, 2) /= n) then
            deallocate(this%f)
            allocate(this%f(m, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%f(m, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_ywork)) then
        if (size(this%m_ywork) /= neqn) then
            deallocate(this%m_ywork)
            allocate(this%m_ywork(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_ywork(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%a)) then
        if (size(this%a, 1) /= n .or. size(this%a, 2) /= n) then
            deallocate(this%a)
            allocate(this%a(n, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%a(n, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%b)) then
        if (size(this%b) /= n) then
            deallocate(this%b)
            allocate(this%b(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%b(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%c)) then
        if (size(this%c) /= n) then
            deallocate(this%c)
            allocate(this%c(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%c(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%e)) then
        if (size(this%e) /= n) then
            deallocate(this%e)
            allocate(this%e(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%e(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! Define the model parameters
    call this%define_model()

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(errmgr, "rkv_alloc_workspace", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine rkv_reset(this)
    !! Resets the integrator to its initial state.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    call this%set_is_first_step(.true.)
end subroutine

! ------------------------------------------------------------------------------
subroutine rkv_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
    fprev, err)
    !! Attempts a single integration step.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
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
        !! An N-element array where the values of the error estimates wil
        !! be written.
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
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.

    ! Local Variables
    integer(int32) :: i, j, n, neqn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = this%get_stage_count()
    neqn = size(y)

    ! Ensure the workspace is allocated
    if (.not.allocated(this%f) .or. .not.allocated(this%m_ywork)) then
        call this%initialize(neqn, errmgr)
        if (errmgr%has_error_occurred()) return
    end if

    ! The Butcher tableau is lower triangular as this is an explicit integrator
    if (.not.this%is_fsal() .or. this%get_is_first_step()) then
        ! On FSAL integrators, we only need to make this call on the first step
        ! as the integrator uses the last evaluation from the previous step
        ! as this step.  On non-FSAL integrators we always need to compute an
        ! updated first step.
        call sys%ode(x, y, this%f(:,1))
    end if
    do i = 2, n
        this%m_ywork = 0.0d0
        do j = 1, i - 1 ! only reference the sub-diagonal components
            this%m_ywork = this%m_ywork + this%a(i, j) * &
                this%f(:,j)
        end do

        call sys%ode( &
            x + h * this%c(i), &
            y + h * this%m_ywork, &
            this%f(:,i) &  ! output
        )
    end do

    ! Compute the two solution estimates, and the resulting error estimate
    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%b(i) * this%f(:,i)
        else
            this%m_ywork = this%m_ywork + this%b(i) * &
                this%f(:,i)
        end if
    end do
    yn = y + h * this%m_ywork

    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%e(i) * this%f(:,i)
        else
            this%m_ywork = this%m_ywork + this%e(i) * &
                this%f(:,i)
        end if
    end do
    en = h * this%m_ywork
end subroutine

! ------------------------------------------------------------------------------
subroutine rkv_on_successful_step(this, x, xn, y, yn)
    !! Perform necessary actions on completion of a successful step.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in) :: xn
        !! The value of the independent variable at the next step.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution values.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution values at the next step.

    ! Local Variables
    integer(int32) :: n

    ! Set up the interpolation polynomial - TO DO check for dense output first
    call this%set_up_interpolation(x, xn, y, yn, this%f)

    ! Store the last result as the first, if this is FSAL
    if (this%is_fsal()) then
        call this%set_is_first_step(.false.)
        n = this%get_stage_count()
        this%f(:,1) = this%f(:,n)
    end if
end subroutine


! ------------------------------------------------------------------------------
pure function rkv_get_alpha(this) result(rst)
    !! Gets the \( \alpha \) control parameter in the PI controller
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
    !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
    !! factor, and \f$ k \) is the order of the integration method.
    class(rk_variable_integrator), intent(in) :: this
        !! The rk_variable_integrator object.
    real(real64) :: rst
        !! The parameter.
    rst = this%m_alpha
end function

! --------------------
subroutine rkv_set_alpha(this, x)
    !! Sets the \( \alpha \) control parameter in the PI controller
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
    !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
    !! factor, and \f$ k \) is the order of the integration method.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    real(real64) :: x
        !! The parameter.
    this%m_alpha = x
end subroutine

! ------------------------------------------------------------------------------
pure function rkv_get_beta(this) result(rst)
    !! Gets the \( \beta \) control parameter in the PI controller
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
    !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
    !! factor, and \f$ k \) is the order of the integration method.
    class(rk_variable_integrator), intent(in) :: this
        !! The rk_variable_integrator object.
    real(real64) :: rst
        !! The parameter.
    rst = this%m_beta
end function

! --------------------
subroutine rkv_set_beta(this, x)
    !! Sets the \( \beta \) control parameter in the PI controller
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} 
    !! e_n^{\alpha} e_{n-1}^{\beta}\), where \( f \) is a safety 
    !! factor, and \f$ k \) is the order of the integration method.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    real(real64), intent(in) :: x
        !! The parameter.
    this%m_beta = x
end subroutine

! ------------------------------------------------------------------------------
function rkv_next_step(this, hn, en, enm1) result(rst)
    !! Computes the next step size using a PI type controller such that
    !! \( h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
    !! e_{n-1}^{\beta} \) where \( f \) is the safety factor and \( k \)
    !! is the order of the integration method.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    real(real64), intent(in) :: hn
        !! The current step size.
    real(real64), intent(in) :: en
        !! The norm of the error for the current step size.
    real(real64), intent(in) :: enm1
        !! The norm of the error from the previous step size.
    real(real64) :: rst
        !! The new step size.

    ! Process
    rst = pi_controller(this%get_alpha(), this%get_beta(), this%get_order(), &
        hn, en, enm1, this%get_safety_factor(), this%get_max_step_size())
end function


! ------------------------------------------------------------------------------
pure function rkv_get_is_first_step(this) result(rst)
    !! Gets a value determining if the integrator is set up to take its first
    !! integration step.
    class(rk_variable_integrator), intent(in) :: this
        !! The rk_variable_integrator object.
    logical :: rst
        !! True if the integrator is on its first step; else, false.
    rst = this%m_firstStep
end function

! --------------------
subroutine rkv_set_is_first_step(this, x)
    !! Sets a value determining if the integrator is set up to take its first
    !! integration step.
    class(rk_variable_integrator), intent(inout) :: this
        !! The rk_variable_integrator object.
    logical :: x
        !! True if the integrator is on its first step; else, false.
    this%m_firstStep = x
end subroutine

! ******************************************************************************
! DORMAND-PRINCE 4/5 INTEGRATOR ROUTINES
! ------------------------------------------------------------------------------
subroutine dprk45_define_model(this)
    use diffeq_dprk45_constants
    !! Defines (initializes) the model parameters.
    class(dprk45_integrator), intent(inout) :: this
        !! The dprk45_integrator object.

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%a = 0.0d0
    
    this%a(2,1) = a21

    this%a(3,1) = a31
    this%a(3,2) = a32

    this%a(4,1) = a41
    this%a(4,2) = a42
    this%a(4,3) = a43

    this%a(5,1) = a51
    this%a(5,2) = a52
    this%a(5,3) = a53
    this%a(5,4) = a54

    this%a(6,1) = a61
    this%a(6,2) = a62
    this%a(6,3) = a63
    this%a(6,4) = a64
    this%a(6,5) = a65
    
    this%a(7,1) = a71
    this%a(7,2) = a72
    this%a(7,3) = a73
    this%a(7,4) = a74
    this%a(7,5) = a75
    this%a(7,6) = a76

    ! B
    this%b(1) = a71
    this%b(2) = a72
    this%b(3) = a73
    this%b(4) = a74
    this%b(5) = a75
    this%b(6) = a76
    this%b(7) = 0.0d0

    ! C
    this%c(1) = 0.0d0
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4
    this%c(5) = c5
    this%c(6) = c6
    this%c(7) = c7

    ! E
    this%e(1) = e1
    this%e(2) = e2
    this%e(3) = e3
    this%e(4) = e4
    this%e(5) = e5
    this%e(6) = e6
    this%e(7) = e7

    ! Update definition status
    this%m_modelDefined = .true.
end subroutine

! ------------------------------------------------------------------------------
pure function dprk45_is_fsal(this) result(rst)
    !! Determines if the integrator is an FSAL (first same as last)
    !! integrator.
    class(dprk45_integrator), intent(in) :: this
        !! The dprk45_integrator object.
    logical :: rst
        !! Returns true, as this integrator is a FSAL integrator.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function dprk45_get_stage_count(this) result(rst)
    !! Gets the number of stages used by the integrator.
    class(dprk45_integrator), intent(in) :: this
        !! The dprk45_integrator object.
    integer(int32) :: rst
        !! Returns the number of stages, 7 for this integrator.
        
    rst = 7
end function

! ------------------------------------------------------------------------------
pure function dprk45_get_order(this) result(rst)
    !! Returns the order of the integrator.
    class(dprk45_integrator), intent(in) :: this
        !! The dprk45_integrator object.
    integer(int32) :: rst
        !! Returns the order of the integrator, 5 for this integrator.
    rst = 5
end function

! ------------------------------------------------------------------------------
subroutine dprk45_set_up_interp(this, x, xn, y, yn, k)
    use diffeq_dprk45_constants
    !! Sets up the interpolation polynomial.
    class(dprk45_integrator), intent(inout) :: this
        !! The dprk45_integrator object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in) :: xn
        !! The value of the independent variable at the next step.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution
        !! values.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution values at
        !! the next step.
    real(real64), intent(in), dimension(:,:) :: k
        !! An N-by-M matrix containing the intermediate step function outputs 
        !! where M is the number of stages of the integrator.

    ! Parameters
    integer(int32), parameter :: n = 5

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h, ydiff, bspl

    ! Intialization
    neqn = size(y)
    h = this%get_step_size()

    ! Memory Allocation
    if (allocated(this%m_dprk45work)) then
        if (size(this%m_dprk45work, 1) /= neqn .or. &
            size(this%m_dprk45work, 2) /= n) &
        then
            deallocate(this%m_dprk45work)
            allocate(this%m_dprk45work(neqn, n))
        end if
    else
        allocate(this%m_dprk45work(neqn, n))
    end if

    ! Construct the coefficient arrays
    do i = 1, neqn
        ydiff = yn(i) - y(i)
        bspl = h * k(i,1) - ydiff

        this%m_dprk45work(i,1) = y(i)
        this%m_dprk45work(i,2) = ydiff
        this%m_dprk45work(i,3) = bspl
        this%m_dprk45work(i,4) = ydiff - h * k(i,7) - bspl
        this%m_dprk45work(i,5) = h * (d1 * k(i,1) + d3 * k(i,3) + &
            d4 * k(i,4) + d5 * k(i,5) + d6 * k(i,6) + d7 * k(i,7))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine dprk45_interp(this, xprev, yprev, xnew, x, y, err)
    !! Provides interpolation between integration points allowing for dense 
    !! output.
    class(dprk45_integrator), intent(in) :: this
        !! The dprk45_integrator object.
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
    integer(int32) :: i, neqn
    real(real64) :: h, s, s1
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(this%m_dprk45work, 1)
    h = xnew - xprev
    s = (x - xprev) / h
    s1 = 1.0d0 - s

    ! Input Check
    if (size(y) /= neqn) then
        call report_array_size_error(errmgr, "dprk45_interp", "y", neqn, &
            size(y))
        return
    end if

    ! Process
    y = this%m_dprk45work(:,1) + s * ( &
        this%m_dprk45work(:,2) + s1 * ( &
        this%m_dprk45work(:,3) + s * ( &
        this%m_dprk45work(:,4) + s1 * this%m_dprk45work(:,5) &
    )))

    ! End
    return
end subroutine

! ******************************************************************************
! BOGACKI-SHAMPINE 3/2 INTEGRATOR ROUTINES
! ------------------------------------------------------------------------------
subroutine bsrk32_define_model(this)
    use diffeq_bsrk32_constants
    !! Defines (initializes) the model parameters.
    class(bsrk32_integrator), intent(inout) :: this
        !! The bsrk32_integrator object.

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%a = 0.0d0

    this%a(2,1) = a21

    this%a(3,1) = a31
    this%a(3,2) = a32

    this%a(4,1) = a41
    this%a(4,2) = a42
    this%a(4,3) = a43

    ! B
    this%b(1) = b1
    this%b(2) = b2
    this%b(3) = b3
    this%b(4) = b4

    ! C
    this%c(1) = c1
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4

    ! E
    this%e(1) = b1a - b1
    this%e(2) = b2a - b2
    this%e(3) = b3a - b3
    this%e(4) = b4a - b4
end subroutine

! ------------------------------------------------------------------------------
pure function bsrk32_is_fsal(this) result(rst)
    !! Determines if the integrator is an FSAL (first same as last)
    !! integrator.
    class(bsrk32_integrator), intent(in) :: this
        !! The bsrk32_integrator object.
    logical :: rst
        !! Returns true, as this integrator is a FSAL integrator.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function bsrk32_get_stage_count(this) result(rst)
    !! Gets the number of stages used by the integrator.
    class(bsrk32_integrator), intent(in) :: this
        !! The bsrk32_integrator object.
    integer(int32) :: rst
        !! Returns the number of stages, 4 for this integrator.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function bsrk32_get_order(this) result(rst)
    !! Returns the order of the integrator.
    class(bsrk32_integrator), intent(in) :: this
        !! The bsrk32_integrator object.
    integer(int32) :: rst
        !! Returns the order of the integrator, 3 for this integrator.
    rst = 3
end function

! ------------------------------------------------------------------------------
subroutine bsrk32_set_up_interp(this, x, xn, y, yn, k)
    !! Sets up the interpolation polynomial.
    class(bsrk32_integrator), intent(inout) :: this
        !! The bsrk32_integrator object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in) :: xn
        !! The value of the independent variable at the next step.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution
        !! values.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution values at
        !! the next step.
    real(real64), intent(in), dimension(:,:) :: k
        !! An N-by-M matrix containing the intermediate step function outputs 
        !! where M is the number of stages of the integrator.

    ! Parameters
    integer(int32), parameter :: n = 3

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h

    ! Initialization
    neqn = size(y)
    h = this%get_step_size()

    ! Memory ALlocation
    if (allocated(this%m_bsrk23work)) then
        if (size(this%m_bsrk23work, 1) /= neqn .or. &
            size(this%m_bsrk23work, 2) /= n) &
        then
            deallocate(this%m_bsrk23work)
            allocate(this%m_bsrk23work(neqn, n))
        end if
    else
        allocate(this%m_bsrk23work(neqn, n))
    end if

    ! Construct the coefficient arrays
    this%m_bsrk23work(:,1) = -(y - yn + k(:,1) * h) / h**2
    this%m_bsrk23work(:,2) = k(:,1) - 2.0d0 * x * this%m_bsrk23work(:,1)
    this%m_bsrk23work(:,3) = y - this%m_bsrk23work(:,1) * x**2 - &
        this%m_bsrk23work(:,2) * x
end subroutine

! ------------------------------------------------------------------------------
subroutine bsrk32_interp(this, xprev, yprev, xnew, x, y, err)
    !! Provides interpolation between integration points allowing for dense 
    !! output.
    class(bsrk32_integrator), intent(in) :: this
        !! The bsrk32_integrator object.
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
    integer(int32) :: neqn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(this%m_bsrk23work, 1)

    ! Input Check
    if (size(y) /= neqn) then
        call report_array_size_error(errmgr, "bsrk32_interp", "y", neqn, &
            size(y))
        return
    end if

    ! Process
    y = this%m_bsrk23work(:,1) * x**2 + this%m_bsrk23work(:,2) * x + &
        this%m_bsrk23work(:,3)

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
end module