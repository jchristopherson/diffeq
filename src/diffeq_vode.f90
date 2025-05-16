module diffeq_vode
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    use ferror
    implicit none
    private
    public :: VODE_ADAMS_METHOD
    public :: VODE_BDF_METHOD
    public :: vode
    public :: adams
    public :: bdf

    integer(int32), parameter :: VODE_ADAMS_METHOD = 10
        !! Describes the VODE Adams method solver.
    integer(int32), parameter :: VODE_BDF_METHOD = 21
        !! Describes the VODE BDF method solver.

    type, abstract, extends(ode_integrator) :: vode
        !! This type encpsulats the VODE variable coefficient backward 
        !! difference code or Adams method.
    contains
        procedure, public :: solve => vode_solve
        procedure, public :: get_order => vode_order_inquiry
        procedure(vode_integer_inquiry), deferred, public :: get_method
    end type

    type, extends(vode) :: adams
        !! Defines an Adams method solver implemented by VODE.
    contains
        procedure, public :: get_method => adams_method_inquiry
    end type

    type, extends(vode) :: bdf
        !! Defines a BDF solver implemented by VODE.
    contains
        procedure, public :: get_method => bdf_method_inquiry
    end type

    interface
        pure function vode_integer_inquiry(this) result(rst)
            !! Defines the signature of a function for inquiring about an
            !! integer-valued property from a vode object.
            use iso_fortran_env, only : int32
            import vode
            class(vode), intent(in) :: this
                !! The vode object.
            integer(int32) :: rst
                !! The requested integer value.
        end function
    end interface

    type vode_argument_container
        !! A container that can be used to pass function pointers and user
        !! information to the VODE code.
        class(*), pointer :: args
            !! User defined arguments.
        logical :: uses_args
            !! Set to true if args is used; else, false.
        procedure(ode), pointer, nopass :: fcn
            !! A pointer to the ODE routine.
        procedure(ode_jacobian), pointer, nopass :: jac
            !! A pointer to the ODE jacobian routine.
    end type

    interface
        subroutine DVODE(f, neqn, y, t, tout, itol, rtol, atol, itask, istate, &
            iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
            use iso_fortran_env, only : real64, int32
            
            interface
                subroutine f(neqn_, t_, y_, ydot_, rpar_, ipar_)
                    use iso_fortran_env, only : int32, real64
                    integer(int32), intent(in) :: neqn_
                    real(real64), intent(in) :: t_, y_(neqn_)
                    real(real64), intent(out) :: ydot_(neqn_)
                    real(real64), intent(inout) :: rpar_(*)
                    integer(int32), intent(inout) :: ipar_(*)
                end subroutine

                subroutine jac(neqn_, t_, y_, ml_, mu_, pd_, nrowpd_, &
                    rpar_, ipar_)
                    use iso_fortran_env, only : int32, real64
                    integer(int32), intent(in) :: neqn_, ml_, mu_, nrowpd_
                    real(real64), intent(in) :: t_, y_(neqn_)
                    real(real64), intent(out) :: pd_(nrowpd_,neqn_)
                    real(real64), intent(inout) :: rpar_(*)
                    integer(int32), intent(inout) :: ipar_(*)
                end subroutine
            end interface

            integer(int32), intent(in) :: neqn, itol, itask, lrw, liw, mf, iopt
            integer(int32), intent(inout) :: istate, iwork(*), ipar(*)
            real(real64), intent(in) :: tout, rtol, atol
            real(real64), intent(inout) :: y(neqn), t, rwork(*), rpar(*)

        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
subroutine vode_solve(this, sys, x, iv, args, err)
    !! Solves the supplied system of ODE's.
    class(vode), intent(inout) :: this
        !! The vode object.
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
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.
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
    integer(int32) :: i, ipar(1), itol, itask, istate, iopt, lrw, liw, mf, nx, &
        neqn, flag, maxord, lwm, miter, nsteps, j, stepsTaken, netf, ncfn
    integer(int32), allocatable, dimension(:) :: iwork
    real(real64) :: rtol, atol, t, tout, xmax
    real(real64), allocatable, dimension(:) :: rpar, rwork, y
    type(vode_argument_container) :: container
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nx = size(x)
    xmax = x(nx)
    neqn = size(iv)
    nsteps = this%get_step_limit()

    ! Input Checking
    if (nx < 2) then
        call report_array_size_error(errmgr, "vode_solve", "x", 2, nx)
        return
    end if
    if (.not.sys%get_is_ode_defined()) then
        call report_missing_ode(errmgr, "vode_solve")
        return
    end if

    ! Additional Initialization
    container%fcn => sys%fcn
    if (associated(sys%jacobian)) container%jac => sys%jacobian
    container%uses_args = present(args)
    if (present(args)) container%args => args
    rpar = transfer(container, rpar)
    ipar(1) = size(rpar)
    allocate(y(neqn), source = iv, stat = flag)
    if (flag /= 0) go to 10
    itol = 1
    rtol = this%get_relative_tolerance()
    atol = this%get_absolute_tolerance()
    istate = 1
    iopt = 1
    miter = 1
    if (this%get_method() == VODE_ADAMS_METHOD) then
        mf = 10
    else
        if (associated(sys%jacobian)) then
            mf = 21
        else
            mf = 22
        end if
    end if
    maxord = this%get_order()

    ! Determine the proper task
    if (nx == 2) then
        if (this%get_allow_overshoot()) then
            itask = 2
        else
            itask = 5
        end if
    else
        if (this%get_allow_overshoot()) then
            itask = 1
        else
            itask = 4
        end if
    end if

    ! Workspace Initializations
    lwm = 2 * neqn**2 + 2
    lrw = 20 + neqn * (maxord + 1) + 3 * neqn + lwm
    liw = 30 + neqn
    allocate(iwork(liw), source = 0, stat = flag)
    if (flag /= 0) go to 10
    allocate(rwork(lrw), source = 0.0d0, stat = flag)
    if (flag /= 0) go to 10

    ! Optional Parameter Initializations
    rwork(1) = xmax
    rwork(6) = this%get_maximum_step_size()
    rwork(7) = this%get_minimum_step_size()

    ! Process
    t = x(1)
    if (nx == 2) then
        tout = x(nx)
    else
        tout = x(2)
    end if
    j = 1
    call this%append_to_buffer(t, y, err = errmgr)
    if (errmgr%has_error_occurred()) return
    do i = 1, nsteps
        ! Take the step
        call DVODE(vode_eqn, neqn, y, t, tout, itol, rtol, atol, itask, &
            istate, iopt, rwork, lrw, iwork, liw, vode_jacobian, mf, rpar, &
            ipar)

        ! Check the state of the integrator
        select case (istate)
        case (1)
            ! Nothing was done as t == tout
        case (-1)
            ! To much work (more than mxstep)
            stepsTaken = iwork(11)
            call report_excessive_iterations(errmgr, "vode_solve", &
                stepsTaken, t)
            return
        case (-2)
            ! Tolerance values are too small
            call report_tolerance_too_small_error(errmgr, "vode_solve")
            return
        case (-3)
            ! Illegal input
            call errmgr%report_error("vode_solve", &
                "An invalid argument was passed to DVODE.  See the " // &
                "command line output for more information.", &
                DIFFEQ_INVALID_INPUT_ERROR)
            return
        case (-4)
            ! Repeated error test failures
            netf = iwork(22)
            call report_successive_error_test_failures(err, "vode_solve", netf)
            return
        case (-5)
            ! Failure to converge
            ncfn = iwork(21)
            call report_multiple_convergence_error(errmgr, "vode_solve", ncfn)
            return
        case (-6)
            ! Pure relative error control requested but failed
            call errmgr%report_error("vode_solve", &
                "Pure relative error control requested, but is " // &
                "unsuccessful for this problem.", DIFFEQ_ERROR_TEST_FAILURE)
            return
        end select

        ! Store the results
        call this%append_to_buffer(t, y, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Update TOUT if nx /= 2
        if (nx /= 2) then
            j = j + 1
            if (j >= nx) go to 100
            tout = x(j)
        end if

        ! Are we done?
        if (abs(t) >= abs(xmax)) then
            ! We're done
            go to 100
        end if
    end do

    ! End
100 continue
    return

10  continue
    ! Memory Error Handling
    call report_memory_error(errmgr, "vode_solve", flag)
    return
end subroutine

! ------------------------------------------------------------------------------
pure function vode_order_inquiry(this) result(rst)
    !! Gets the highest order of this integrator.  This integrator is a variable
    !! order integrator; therefore, the exact order is problem dependent.
    class(vode), intent(in) :: this
        !! The vode object.
    integer(int32) :: rst
        !! The highest order of this integrator.

    ! Process
    if (this%get_method() == VODE_ADAMS_METHOD) then
        rst = 12
    else
        rst = 5
    end if
end function

! ------------------------------------------------------------------------------
subroutine vode_eqn(neqn, t, y, ydot, rpar, ipar)
    !! The routine containing the differential equations solved by VODE.
    integer(int32), intent(in) :: neqn
        !! The number of equations.
    real(real64), intent(in) :: t
        !! The current value of the independent variable.
    real(real64), intent(in) :: y(neqn)
        !! The current state vector.
    real(real64), intent(out) :: ydot(neqn)
        !! The output derivative vector.
    real(real64), intent(inout) :: rpar(*)
        !! Real-valued parameter array for communication with the calling code.
        !! This is the array used to transfer data.  The ipar array is used
        !! to transfer information regarding size.
    integer(int32), intent(inout) :: ipar(*)
        !! Integer-valued parameter array for communication with the calling
        !! code.  Only the first element of this array is used; it contains
        !! the size of the stored array rpar.

    ! Local Variables
    type(vode_argument_container) :: args

    ! Extract information from the user
    args = transfer(rpar(1:ipar(1)), args)

    ! Evaluate the routine
    if (args%uses_args) then
        call args%fcn(t, y(1:neqn), ydot(1:neqn), args%args)
    else
        call args%fcn(t, y(1:neqn), ydot(1:neqn))
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine vode_jacobian(neqn, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    !! The routine containing the Jacobian calculation routine.
    integer(int32), intent(in) :: neqn
        !! The number of equations.
    real(real64), intent(in) :: t
        !! The current value of the independent variable.
    real(real64), intent(in) :: y(neqn)
        !! The current state vector.
    integer(int32), intent(in) :: ml
        !! The lower bandwidth of a banded Jacobian.
    integer(int32), intent(in) :: mu
        !! The upper bandwidth of a banded Jacobian.
    integer(int32), intent(in) :: nrowpd
        !! The leading dimension of the Jacobian matrix.
    real(real64), intent(out) :: pd(nrowpd,neqn)
        !! The Jacobian matrix.

    real(real64), intent(inout) :: rpar(*)
        !! Real-valued parameter array for communication with the calling code.
        !! This is the array used to transfer data.  The ipar array is used
        !! to transfer information regarding size.
    integer(int32), intent(inout) :: ipar(*)
        !! Integer-valued parameter array for communication with the calling
        !! code.  Only the first element of this array is used; it contains
        !! the size of the stored array rpar.

    ! Local Variables
    type(vode_argument_container) :: args

    ! Extract information from the user
    args = transfer(rpar(1:ipar(1)), args)

    ! Evaluate the routine
    if (args%uses_args) then
        call args%jac(t, y(1:neqn), pd(1:neqn,1:neqn), &
            args%args)
    else
        call args%jac(t, y(1:neqn), pd(1:neqn,1:neqn))
    end if
end subroutine

! ******************************************************************************
! ADAMS
! ------------------------------------------------------------------------------
pure function adams_method_inquiry(this) result(rst)
    !! Gets the method identifier for this integrator.
    class(adams), intent(in) :: this
        !! The adams object.
    integer(int32) :: rst
        !! The method identifier.

    rst = VODE_ADAMS_METHOD
end function

! ******************************************************************************
! BDF
! ------------------------------------------------------------------------------
pure function bdf_method_inquiry(this) result(rst)
    !! Gets the method identifier for this integrator.
    class(bdf), intent(in) :: this
        !! The bdf object.
    integer(int32) :: rst
        !! The method identifier.

    rst = VODE_BDF_METHOD
end function

! ------------------------------------------------------------------------------
end module