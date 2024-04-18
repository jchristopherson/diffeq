module diffeq_base
    !! A collection of base types for the DIFFEQ library.
    use iso_fortran_env
    use diffeq_errors
    use ferror
    implicit none
    private
    public :: ode
    public :: ode_jacobian
    public :: ode_mass_matrix
    public :: ode_container
    public :: ode_integrator
    public :: ode_solver
    public :: ode_integer_inquiry

! ------------------------------------------------------------------------------
    interface
        subroutine ode(x, y, dydx)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        subroutine ode_jacobian(x, y, jac)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        subroutine ode_mass_matrix(x, y, m)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: m
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type ode_container
        !! A container for the routine containing the ODEs to integrate.
        logical, private :: m_massDependent = .true.
            ! A value determining if the mass matrix is state dependent such 
            !! that it must be recomputed at each step. 
        real(real64), private, allocatable, dimension(:) :: m_jwork
            ! Jacobian calculation workspace array.
        real(real64), private :: m_fdStep = sqrt(epsilon(1.0d0))
            ! Finite difference step size.
        procedure(ode), pointer, public, nopass :: fcn => null()
            !! A pointer to the routine containing the ODEs to integrate.
        procedure(ode_jacobian), pointer, public, nopass :: &
            jacobian => null()
            !! A pointer to the routine containing the analytical Jacobian.
            !! If supplied, this routine is utilized; however, if null, a finite
            !! difference approximation is utilized.
        procedure(ode_mass_matrix), pointer, public, nopass :: &
            mass_matrix => null()
            !! A pointer to the routine containing the mass matrix for the
            !! system.  If set to null (the default), an identity mass matrix 
            !! will be assumed.
    contains
        procedure, private :: allocate_workspace => oc_alloc_workspace
        procedure, public :: get_finite_difference_step => oc_get_fd_step
        procedure, public :: set_finite_difference_step => oc_set_fd_step
        procedure, public :: get_is_mass_matrix_dependent => &
            oc_get_is_mass_dependent
        procedure, public :: set_is_mass_matrix_dependent => &
            oc_set_is_mass_dependent
        procedure, public :: compute_jacobian => oc_jacobian
        procedure, public :: ode => oc_ode_fcn
        procedure, public :: get_is_ode_defined => oc_get_is_ode_defined
    end type

! ------------------------------------------------------------------------------
    type, abstract :: ode_integrator
        !! The most basic ODE integrator object capable of integrating
        !! systems of ODE's.
    contains
        procedure(ode_solver), public, pass, deferred :: solve
            !! Solves the supplied system of ODE's.
        procedure(ode_integer_inquiry), public, pass, deferred :: get_order
            !! Returns the order of the integrator.
    end type

    interface
        function ode_solver(this, sys, x, iv, err) result(rst)
            !! Solves the supplied system of ODE's.
            use iso_fortran_env
            use ferror
            import ode_integrator
            import ode_container
            class(ode_integrator), intent(inout) :: this
                !! The ode_integrator object.
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
        end function

        pure function ode_integer_inquiry(this) result(rst)
            !! Returns an integer value from the ode_integrator object.
            use iso_fortran_env
            import ode_integrator
            class(ode_integrator), intent(in) :: this
                !! The ode_integrator object.
            integer(int32) :: rst
                !! The requested value.
        end function
    end interface


    
contains
! ******************************************************************************
! ODE_CONTAINER ROUTINES
! ------------------------------------------------------------------------------
pure function oc_get_fd_step(this) result(rst)
    !! Gets the size of the step to use for the finite difference
    !! calculations used to estimate the Jacobian.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    real(real64) :: rst
        !! The step size.
    rst = this%m_fdStep
end function

! --------------------
subroutine oc_set_fd_step(this, x)
    !! Sets the size of the step to use for the finite difference
    !! calculations used to estimate the Jacobian.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    real(real64), intent(in) :: x
        !! The step size.
    this%m_fdStep = x
end subroutine

! ------------------------------------------------------------------------------
pure function oc_get_is_mass_dependent(this) result(rst)
    !! Gets a value determining if the mass matrix is state-dependent
    !! such that it requires updating at every integration step.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    logical :: rst
        !! True if the mass matrix is state-dependent such that it 
        !! requires updating at each integration step; else, false if the
        !! mass matrix is not state-dependent and can be treated as constant
        !! for all integration steps.
    rst = this%m_massDependent
end function

! --------------------
subroutine oc_set_is_mass_dependent(this, x)
    !! Sets a value determining if the mass matrix is state-dependent
    !! such that it requires updating at every integration step.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    logical :: x
        !! True if the mass matrix is state-dependent such that it 
        !! requires updating at each integration step; else, false if the
        !! mass matrix is not state-dependent and can be treated as constant
        !! for all integration steps.
    this%m_massDependent = x
end subroutine

! ------------------------------------------------------------------------------
subroutine oc_jacobian(this, x, y, jac, err)
    !! Computes the Jacobian matrix for the system of ODEs.  If
    !! a routine is provided with an analytical Jacobian, the supplied
    !! routine is utilized; else, the Jacobian is estimated via a forward
    !! difference approximation.
    class(ode_container), intent(inout) :: this
        !! The ode_container object.
    real(real64), intent(in) :: x
        !! The current independent variable value.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current dependent
        !! variable values.
    real(real64), intent(out), dimension(:,:) :: jac
        !! An N-by-N matrix where the Jacobian will be written.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling. Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined,
        !!      and the calculation is being performed by finite differences.
        !! 
        !!  - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if jac is not N-by-N.

    ! Local Variables
    integer(int32) :: i, ndof
    real(real64) :: h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    ndof = size(y)
    h = this%get_finite_difference_step()

    ! Input Checking
    if (size(jac, 1) /= ndof .or. size(jac, 2) /= ndof) then
        call report_matrix_size_error(errmgr, "oc_jacobian", "jac", &
            ndof, ndof, size(jac, 1), size(jac, 2))
        return
    end if

    ! Use a user-defined routine, and then be done
    if (associated(this%jacobian)) then
        call this%jacobian(x, y, jac)
        return
    end if

    ! Allocate workspace.  No action is taken if the proper workspace is
    ! already allocated.
    call this%allocate_workspace(ndof, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Finite Difference Approximation
    ! J(i,j) = df(i) / dy(j)
    this%m_jwork(1:ndof) = y
    call this%ode(x, y, this%m_jwork(ndof+1:))
    do i = 1, ndof
        this%m_jwork(i) = this%m_jwork(i) + h
        call this%ode(x, this%m_jwork(1:ndof), jac(:,i))
        jac(:,i) = (jac(:,i) - this%m_jwork(ndof+1:)) / h
        this%m_jwork(i) = y(i)
    end do

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine oc_alloc_workspace(this, ndof, err)
    ! Use to allocate internal workspaces.  This routine only takes action
    ! if the workspace array(s) are not sized properly for the application.
    class(ode_container), intent(inout) :: this
        ! The ode_container object.
    integer(int32), intent(in) :: ndof
        ! The number of degrees of freedom.
    class(errors), intent(inout) :: err
        ! The error handling object.

    ! Local Variables
    integer(int32) :: flag

    ! Jacobian Workspace Allocation
    if (allocated(this%m_jwork)) then
        if (size(this%m_jwork) /= 2 * ndof) then
            deallocate(this%m_jwork)
            allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
            if (flag /= 0) then
                call report_memory_error(err, "oc_alloc_workspace", flag)
                return
            end if
        end if
    else
        allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
        if (flag /= 0) then
            call report_memory_error(err, "oc_alloc_workspace", flag)
            return
        end if
    end if

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine oc_ode_fcn(this, x, y, dydx)
    !! Evaluates the ODEs by evaluating the routine defined by
    !! fcn.  This routine may also be overidden to provide custom
    !! functionallity.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current values of the
        !! dependent variables.
    real(real64), intent(out), dimension(:) :: dydx
        !! An N-element array where the output of each of the
        !! N ODEs will be written.

    ! Process
    if (associated(this%fcn)) then
        call this%fcn(x, y, dydx)
    end if
end subroutine

! ------------------------------------------------------------------------------
pure function oc_get_is_ode_defined(this) result(rst)
    !! Gets a logical value determining if the ODE routine has been defined.
    class(ode_container), intent(in) :: this
        !! The ode_container object.
    logical :: rst
        !! True if the ODE routine has been defined; else, false.

    rst = associated(this%fcn)
end function

! ------------------------------------------------------------------------------
end module