module diffeq_fixed_step
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    implicit none
    private
    public :: fixed_step_integrator
    public :: ode_fixed_step

! ------------------------------------------------------------------------------
    type, abstract, extends(ode_integrator) :: fixed_step_integrator
        !! Defines a fixed-step integrator.
    contains
        procedure, public :: solve => fsi_solver
            !! Solves the supplied system of ODE's.
        procedure(ode_fixed_step), public, pass, deferred :: step
            !! Takes a single integration step.
    end type

    ! diffeq_fs_integrator.f90
    interface
        subroutine ode_fixed_step(this, sys, h, x, y, yn, xprev, yprev, fprev, &
            err)
            !! Takes a single fixed-size integration step.
            use iso_fortran_env 
            use ferror   
            import fixed_step_integrator
            import ode_container
            class(fixed_step_integrator), intent(inout) :: this
                !! The fixed_step_integrator object.
            class(ode_container), intent(inout) :: sys
                !! The ode_container object containing the ODE's to integrate.
            real(real64), intent(in) :: h
                !! The size of the step to take.
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current values of the
                !! dependent variables.
            real(real64), intent(out), dimension(:) :: yn
                !! An N-element array where the values of the dependent 
                !! variables at x + h will be written.
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
                !! implementation of the errors class is used internally to 
                !! provide error handling.
        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
function fsi_solver(this, sys, x, iv, err) result(rst)
    !! Solves the supplied system of ODE's.
    class(fixed_step_integrator), intent(inout) :: this
        !! The fixed_step_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array, of at least 2 values, defining at a minimum
        !! the starting and ending values of the independent variable 
        !! integration range.  If more than two values are specified, 
        !! the integration results will be returned at the supplied 
        !! values.
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
        !! - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
        !!      memory allocation issue.
        !!
        !! - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is 
        !!      defined.
        !!
        !! - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if there are less than 
        !!      2 values given in the independent variable array x.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An M-by-N matrix where M is the number of solution points, 
        !! and N is the number of ODEs plus 1.  The first column 
        !! contains the values of the independent variable at which the 
        !! results were computed.  The remaining columns contain the 
        !! integration results for each ODE.

    ! Local Variables
    integer(int32) :: i, j, npts, neqn, flag
    real(real64) :: h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    neqn = size(iv)

    ! Input Checking
    if (npts < 2) then
        call report_min_array_size_not_met(errmgr, "fsi_solver", "x", 2, npts)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "fsi_solver")
        return
    end if

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "fsi_solver", flag)
        return
    end if

    ! Process
    rst(1,1) = x(1)
    rst(1,2:) = iv
    do i = 2, npts
        ! Compute the solution
        j = i - 1
        h = x(i) - x(j)
        rst(i,1) = x(i)
        call this%step(sys, h, rst(j,1), rst(j,2:), rst(i,2:), err = errmgr)
        if (errmgr%has_error_occurred()) return
    end do

    ! End
    return
end function

! ------------------------------------------------------------------------------
end module