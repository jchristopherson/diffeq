module diffeq_multistep_fixed
    use iso_fortran_env
    use diffeq_fixed_step
    use diffeq_errors
    use diffeq_rk4
    use diffeq_base
    implicit none
    private
    public :: fixed_multistep_integrator

    type, abstract, extends(fixed_step_integrator) :: fixed_multistep_integrator
        !! Defines a fixed step-size, multi-step integrator.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            ! An NEQN-by-ORDER storage matrix for ODE outputs.
    contains
        procedure, private :: allocate_workspace => fms_alloc_workspace
            !! Allocates workspace variables.
        procedure, public :: solve => fms_solver
            !! Solves the supplied system of ODEs.
    end type

contains
! ------------------------------------------------------------------------------
subroutine fms_alloc_workspace(this, neqn, err)
    !! Allocates workspace variables.
    class(fixed_multistep_integrator), intent(inout) :: this
        !! The fixed_multistep_integrator object.
    integer(int32), intent(in) :: neqn
        !! The number of equations to integrate.
    class(errors), intent(inout) :: err
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
    integer(int32) :: n, flag

    ! Process
    n = this%get_order()
    if (allocated(this%m_buffer)) then
        if (size(this%m_buffer, 1) /= neqn .or. &
            size(this%m_buffer, 2) /= n) &
        then
            deallocate(this%m_buffer)
            allocate(this%m_buffer(neqn, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_buffer(neqn, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(err, "fms_alloc_workspace", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
function fms_solver(this, sys, x, iv, err) result(rst)
    !! Solves the supplied system of ODEs.
    class(fixed_multistep_integrator), intent(inout) :: this
        !! The fixed_multistep_integrator object.
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
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a 
        !!      memory allocation issue.
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is 
        !!      defined.
        !!
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if there are less than 
        !!      2 values given in the independent variable array x.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An M-by-N matrix where M is the number of solution points, 
        !! and N is the number of ODEs plus 1.  The first column 
        !! contains the values of the independent variable at which the 
        !! results were computed.  The remaining columns contain the 
        !! integration results for each ODE.

    ! Local Variables
    integer(int32) ::i, j, j1, j2, npts, neqn, order, mn, flag
    real(real64) :: h
    type(rk4_fixed_integrator) :: starter
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    neqn = size(iv)
    order = this%get_order()
    mn = min(order, npts)

    ! Input Checking
    if (npts < 2) then
        call report_min_array_size_not_met(errmgr, "fms_solver", "x", 2, npts)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "fms_solver")
        return
    end if

    ! Memory Allocation
    allocate(rst(npts, neqn + 1), stat = flag)
    if (flag /= 0) then
        call report_memory_error(errmgr, "fms_solver", flag)
        return
    end if

    call this%allocate_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Use a 4th order Runge-Kutta integrator to step into the problem
    rst(1,1) = x(1)
    rst(1,2:) = iv
    call sys%ode(x(1), iv, this%m_buffer(:,order))
    j2 = order - 1
    do i = 2, mn
        j = i - 1
        h = x(j) - x(i)
        call starter%step(sys, h, x(j), rst(j,2:), rst(i,2:))
        call sys%ode(x(i), rst(i,2:), this%m_buffer(:,j2))
        j2 = j2 - 1
    end do

    ! Finish the problem with the multi-step method
    j1 = 1
    j2 = order
    do i = order + 1, npts
        h = x(i) - x(j2)
        rst(i,1) = x(i)
        call this%step(sys, h, rst(j2,1), rst(j2,2:), rst(i,2:), &
            rst(j1:j2,1), rst(j1:j2,2:), this%m_buffer, err = errmgr)
        if (errmgr%has_error_occurred()) return
        j1 = j1 + 1
        j2 = j2 + 1
    end do

    ! End
    return
end function

! ------------------------------------------------------------------------------
end module