module diffeq_variable_singlestep
    use iso_fortran_env
    use diffeq_variable_step
    use diffeq_errors
    use diffeq_base
    implicit none
    private
    public :: variable_singlestep_integrator

    !> @brief Defines a variable-step, single-stage integrator.
    type, abstract, extends(variable_step_integrator) :: &
        variable_singlestep_integrator
    contains
        procedure, public :: solve => vssi_solve
            !! Solves the supplied system of ODEs.
        procedure, private :: solve_driver => vssi_solve_driver
        procedure, private :: dense_solve_driver => vssi_dense_solve_driver
    end type

contains
! ------------------------------------------------------------------------------
function vssi_solve(this, sys, x, iv, err) result(rst)
    !! Solves the supplied system of ODEs.
    class(variable_singlestep_integrator), intent(inout) :: this
        !! The variable_singlestep_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array, of at least 2 values, defining at a minimum
        !! the starting and ending values of the independent variable 
        !! integration range.  If more than two values are specified, the
        !! integration results will be returned at the supplied values.
    real(real64), intent(in), dimension(:) :: iv
        !! An array containing the initial values for each ODE.
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
        !!
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if @p x has less than 2 elements.
        !!
        !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size becomes
        !!      too small.
        !!
        !!  - DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR: Occurs if the iteration
        !!      count is exceeded for a single step.
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE routine is defined.
        !!
        !!  - DIFFEQ_INVALID_INPUT_ERROR: Occurs if max(x) - min(x) = 0.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An M-by-N matrix where M is the number of solution points, 
        !! and N is the number of ODEs plus 1.  The first column contains
        !! the values of the independent variable at which the results were
        !! computed.  The remaining columns contain the integration results
        !! for each ODE.

    ! Local Variables
    integer(int32) :: nx
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nx = size(x)

    ! Input Checking
    if (nx < 2) then
        call report_min_array_size_not_met(errmgr, "vssi_solve", "x", 2, nx)
        return
    end if
    if (abs(maxval(x) - minval(x)) < epsilon(1.0d0)) then
        call errmgr%report_error("vssi_solve", "MAX(X) - MIN(X) is zero.", &
            DIFFEQ_INVALID_INPUT_ERROR)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "vssi_solve")
        return
    end if

    ! Process
    if (nx == 2) then
        rst = this%solve_driver(sys, x, iv, errmgr)
    else
        rst = this%dense_solve_driver(sys, x, iv, errmgr)
    end if
    if (errmgr%has_error_occurred()) return

    ! End
    return
end function

! ------------------------------------------------------------------------------
function vssi_solve_driver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(variable_singlestep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout) :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Parameters
    real(real64), parameter :: sml = 1.0d2 * epsilon(1.0d0)

    ! Local Variables
    integer(int32) :: i, neqn, flag, order
    real(real64) :: xn, xmax
    real(real64), allocatable, dimension(:) :: yn, yn1

    ! Initialization
    xmax = x(2)
    xn = x(1)
    neqn = size(iv)
    order = this%get_order()

    ! Memory Allocation
    allocate(yn(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(yn1(neqn), stat = flag)
    if (flag /= 0) then
        call report_memory_error(err, "vssi_solve_driver", flag)
        return
    end if

    ! Store the initial conditions
    call this%buffer_results(x(1), iv, err)
    if (err%has_error_occurred()) return

    ! Provide an initial step size estimate
    call sys%ode(xn, iv, yn1)
    call this%set_next_step_size( &
        this%estimate_first_step_size(xn, xmax, iv, yn1) &
    )

    ! Cycle until complete
    i = 0
    do
        ! Take the step
        call this%step(sys, xn, xmax, yn, yn1, err = err)
        if (err%has_error_occurred()) return

        ! Update xn and yn
        xn = xn + this%get_step_size()
        yn = yn1

        ! Store the solution
        call this%buffer_results(xn, yn, err)
        if (err%has_error_occurred()) return

        ! Break if |XN| > |XMAX|
        if (abs(xn) >= abs(xmax) .or. abs(this%get_step_size()) < sml) exit

        ! Iteration Counter
        i = i + 1
        if (i > this%get_max_integration_step_count()) then
            call report_excessive_integration_steps(err, &
                "vssi_solve_driver", i, xn)
            return
        end if
    end do

    ! Output the results
    rst = this%get_buffer_contents()

    ! End
    return
end function

! ------------------------------------------------------------------------------
function vssi_dense_solve_driver(this, sys, x, iv, err) result(rst)
    ! Arguments
    class(variable_singlestep_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(errors), intent(inout) :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Parameters
    real(real64), parameter :: sml = 1.0d2 * epsilon(1.0d0)

    ! Local Variables
    integer(int32) :: i, j, npts, neqn, flag
    real(real64) :: xn, xmax, xn1, xi
    real(real64), allocatable, dimension(:) :: yn, yn1
    
    ! Initialization
    neqn = size(iv)
    npts = size(x)
    xn = x(1)
    xmax = x(npts)

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag == 0) allocate(yn(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(yn1(neqn), stat = flag)
    if (flag /= 0) then
        call report_memory_error(err, "vssi_solve_driver", flag)
        return
    end if

    ! Store the initial conditions
    rst(1,1) = x(1)
    rst(1,2:) = iv

    ! Cycle until complete
    i = 0
    j = 1
    xi = x(2)
    outer: do
        ! Take a step
        call this%step(sys, xn, xmax, yn, yn1, err = err)
        if (err%has_error_occurred()) return
        xn1 = xn + this%get_step_size()

        ! Interpolate as needed to achieve any intermediary solution points
        do while (abs(xi) <= abs(xn1))
            j = j + 1
            call this%interpolate(xn, yn, xn1, xi, rst(j,2:), err)
            if (err%has_error_occurred()) return
            rst(j,1) = xi
            if (j >= npts) exit outer
            xi = x(j + 1)
        end do

        ! Update xn and yn
        xn = xn1
        yn = yn1

        ! Break if |XN| > |XMAX|
        if (abs(xn) >= abs(xmax) .or. abs(this%get_step_size()) < sml) exit

        ! Iteration Counter
        i = i + 1
        if (i > max(this%get_max_integration_step_count(), npts)) then
            call report_excessive_integration_steps(err, &
                "vssi_solve_driver", i, xn)
            return
        end if
    end do outer

    ! End
    return
end function

! ------------------------------------------------------------------------------
end module