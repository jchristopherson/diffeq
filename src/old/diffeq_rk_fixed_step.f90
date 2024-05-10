module diffeq_rk_fixed_step
    !! This module contains the structure for fixed-step Runge-Kutta 
    !! integrators.
    use iso_fortran_env
    use diffeq_errors
    use diffeq_base
    use diffeq_fixed_step
    implicit none
    private
    public :: rk_fixed_integrator
    public :: rkf_get_array_parameter
    public :: rkf_get_matrix_parameter

    type, abstract, extends(fixed_step_integrator) :: rk_fixed_integrator
        !! Defines an explicit, Runge-Kutta fixed-step integrator.
        real(real64), private, allocatable, dimension(:,:) :: m_work
            ! Workspace matrix.
    contains
        ! Use to allocate internal workspaces.  This routine only takes action
        ! if the workspace array(s) are not sized properly for the application.
        procedure, private :: allocate_workspace => rkf_alloc_workspace
        procedure, public :: step => rkf_step
            !! Takes a single Runge-Kutta integration step.
        procedure(rkf_get_matrix_parameter), deferred, public :: &
            get_method_factor
            !! Gets the requested method factor from the Butcher tableau.
        procedure(rkf_get_array_parameter), deferred, public :: &
            get_quadrature_weight
            !! Gets the requested quadrature weight from the Butcher tableau.
        procedure(rkf_get_array_parameter), deferred, public :: &
            get_position_factor
            !! Gets the requested position factor from the Butcher tableau.
    end type

    interface
        pure function rkf_get_matrix_parameter(this, i, j) result(rst)
            !! Retrieves a parameter from a matrix stored by the
            !! rk_fixed_integrator object.
            use iso_fortran_env
            import rk_fixed_integrator
            class(rk_fixed_integrator), intent(in) :: this
                !! The rk_fixed_integrator object.
            integer(int32), intent(in) :: i
                !! The row index of the matrix parameter to retrieve.
            integer(int32), intent(in) :: j
                !! The column index of the matrix parameter to retrieve.
            real(real64) :: rst
                !! The requested parameter.
        end function

        pure function rkf_get_array_parameter(this, i) result(rst)
            !! Retrieves a parameter from an array stored by the 
            !! rk_fixed_integrator object.
            use iso_fortran_env
            import rk_fixed_integrator
            class(rk_fixed_integrator), intent(in) :: this
                !! The rk_fixed_integrator object.
            integer(int32), intent(in) :: i
                !! The index of the array parameter to retrieve.
            real(real64) :: rst
                !! The requested parameter
        end function
    end interface


contains
! ------------------------------------------------------------------------------
! n = method order
! neqn = # of ODE's to solve
subroutine rkf_alloc_workspace(this, n, neqn, err)
    ! Arguments
    class(rk_fixed_integrator), intent(inout) :: this
    integer(int32), intent(in) :: n, neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: flag, mw, nw

    ! Process
    mw = neqn
    nw = n + 1  ! +1 for 1 additional NEQN-element workspace array
    if (allocated(this%m_work)) then
        if (size(this%m_work, 1) /= mw .or. size(this%m_work, 2) /= nw) then
            deallocate(this%m_work)
            allocate(this%m_work(mw, nw), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_work(mw, nw), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    call report_memory_error(err, "rkf_alloc_workspace", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine rkf_step(this, sys, h, x, y, yn, xprev, yprev, fprev, err)
    !! Takes a single Runge-Kutta integration step.
    class(rk_fixed_integrator), intent(inout) :: this
        !! The rk_fixed_integrator object.
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
        !! provide error handling.  Possible error and warning messages that
        !! may be encountered are as follows.
        !!
        !! - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if yn is not sized appropriately.

    ! Local Variables
    integer(int32) :: i, j, n, neqn, m, n1
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    n = this%get_order()
    n1 = n + 1  ! index of additional workspace array #1

    ! Input Checking
    if (size(yn) /= neqn) then
        call report_array_size_error(errmgr, "rkf_step", "yn", neqn, size(yn))
        return
    end if

    ! Allocate the workspace
    call this%allocate_workspace(n, neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! As this is an explicit routine, the Butcher tableau is lower triangular.
    call sys%ode(x, y, this%m_work(:,1))
    do i = 2, n
        this%m_work(:,n1) = 0.0d0
        do j = 1, i - 1 ! only reference the sub-diagonal components
            this%m_work(:,n1) = this%m_work(:,n1) + &
                this%get_method_factor(i,j) * this%m_work(:,j)
        end do

        call sys%ode( &
            x + h * this%get_position_factor(i), &
            y + h * this%m_work(:,n1), &
            this%m_work(:,i) &   ! output
        )
    end do

    ! Compute the next solution estimate
    do i = 1, n
        if (i == 1) then
            this%m_work(:,n1) = this%get_quadrature_weight(i) * this%m_work(:,i)
        else
            this%m_work(:,n1) = this%m_work(:,n1) + &
                this%get_quadrature_weight(i) * this%m_work(:,i)
        end if
    end do
    yn = y + h * this%m_work(:,n1)

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
end module