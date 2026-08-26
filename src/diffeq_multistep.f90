module diffeq_multistep
    !! Variable-order multistep ODE solvers backed by DVODE.
    !!
    !! VODE uses variable-step formulas of the general form
    !! \[
    !! \sum_{j=0}^{q} \alpha_j y_{n-j}
    !!   = h_n \beta f(x_n,y_n),
    !! \]
    !! where the coefficient set and order \(q\) change with the selected
    !! method.  Adams formulas are intended primarily for non-stiff systems;
    !! backward differentiation formulas are intended primarily for stiff
    !! systems.  The modern library types manage DVODE's work arrays,
    !! tolerances, callbacks, interpolation, and solution buffering.
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    implicit none
    private
    public :: DIFFEQ_ADAMS_METHOD
    public :: DIFFEQ_BDF_METHOD
    public :: multistep_integrator
    public :: adams
    public :: bdf

    integer(int32), parameter :: DIFFEQ_ADAMS_METHOD = 10
        !! Describes the VODE Adams method solver.
    integer(int32), parameter :: DIFFEQ_BDF_METHOD = 21
        !! Describes the VODE BDF method solver.

    type, abstract, extends(ode_integrator) :: multistep_integrator
        !! Abstract base for the DVODE-backed Adams and BDF solver variants.
        !!
        !! VODE controls a weighted local error against the requested absolute
        !! and relative tolerances, using a component scale such as
        !! \( \mathrm{sc}_i = \mathrm{atol}_i + \mathrm{rtol}_i |y_i| \).
    contains
        procedure, public :: solve => vode_solve
        procedure, public :: get_order => vode_order_inquiry
        procedure(vode_integer_inquiry), deferred, public :: get_method
    end type

    type, extends(multistep_integrator) :: adams
        !! Variable-order Adams solver implemented by VODE.
        !!
        !! Adams predictor--corrector formulas use past derivative information
        !! and are efficient for smooth, non-stiff systems.  They are usually
        !! the preferred VODE method when stiffness is not expected.
    contains
        procedure, public :: get_method => adams_method_inquiry
    end type

    type, extends(multistep_integrator) :: bdf
        !! Variable-order backward differentiation formula solver implemented
        !! by VODE.
        !!
        !! BDF formulas are implicit and use past solution values.  They are
        !! designed for stiff systems and can remain stable at step sizes that
        !! would make an explicit method impractical.
    contains
        procedure, public :: get_method => bdf_method_inquiry
    end type

    interface
        pure function vode_integer_inquiry(this) result(rst)
            !! Defines the signature of a function for inquiring about an
                !! integer-valued property from a multistep integrator.
            use iso_fortran_env, only : int32
            import multistep_integrator
            class(multistep_integrator), intent(in) :: this
                !! The multistep integrator.
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
        procedure(ode_mass_matrix), pointer, nopass :: mass => null()
            !! A pointer to the optional mass matrix routine.
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

        subroutine DVINDY(tout, k, yh, ldyh, y, iflag)
            use iso_fortran_env, only : real64, int32
            integer(int32), intent(in) :: k, ldyh
            integer(int32), intent(out) :: iflag
            real(real64), intent(in) :: tout, yh(ldyh,*)
            real(real64), intent(out) :: y(*)
        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
subroutine vode_solve(this, sys, x, iv, args)
    !! Solves the supplied system of ODE's.
    !!
    !! DVODE controls the variable step size and method order internally.
    !! When `x` contains only the initial and final coordinates, every
    !! successful internal DVODE step is appended to the solution buffer.
    !! When additional coordinates are supplied, DVODE is advanced to those
    !! coordinates and its dense-output interpolation is used as requested.
    class(multistep_integrator), intent(inout) :: this
        !! The multistep integrator.
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

    ! Local Variables
    integer(int32) :: i, ipar(1), itol, itask, istate, iopt, lrw, liw, mf, nx, &
        neqn, maxord, lwm, miter, nsteps, j, stepsTaken, netf, ncfn
    integer(int32), allocatable, dimension(:) :: iwork
    real(real64) :: rtol, atol, t, tout, xmax
    real(real64), allocatable, dimension(:) :: rpar, rwork, y
    type(vode_argument_container) :: container
    
    ! Initialization
    nx = size(x)
    xmax = x(nx)
    neqn = size(iv)
    nsteps = this%get_step_limit()

    ! Input Checking
    if (nx < 2) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (.not.sys%get_is_ode_defined()) error stop DIFFEQ_MISSING_ARGUMENT_ERROR

    ! Package the library callbacks for DVODE's legacy calling convention.
    container%fcn => sys%fcn
    if (associated(sys%jacobian)) container%jac => sys%jacobian
    if (associated(sys%mass_matrix)) container%mass => sys%mass_matrix
    container%uses_args = present(args)
    if (present(args)) container%args => args
    rpar = transfer(container, rpar)
    ipar(1) = size(rpar)
    allocate(y(neqn), source = iv)
    itol = 1
    rtol = this%get_relative_tolerance()
    atol = this%get_absolute_tolerance()
    istate = 1
    iopt = 1
    miter = 1
    if (this%get_method() == DIFFEQ_ADAMS_METHOD) then
        mf = 10
    else
        if (associated(sys%jacobian)) then
            mf = 21
        else
            mf = 22
        end if
    end if
    maxord = this%get_order()

    ! Use one-step output for a two-point request so accepted internal steps
    ! remain visible; DVINDY restores the requested endpoint after overshoot.
    if (nx == 2) then
        itask = 2
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
    allocate(iwork(liw), source = 0)
    allocate(rwork(lrw), source = 0.0d0)

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
    call this%append_to_buffer(t, y)
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
            error stop DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR
        case (-2)
            ! Tolerance values are too small
            error stop DIFFEQ_TOLERANCE_TOO_SMALL
        case (-3)
            ! Illegal input
            error stop DIFFEQ_INVALID_INPUT_ERROR
        case (-4)
            ! Repeated error test failures
            netf = iwork(22)
            error stop DIFFEQ_ERROR_TEST_FAILURE
        case (-5)
            ! Failure to converge
            ncfn = iwork(21)
            error stop DIFFEQ_CONVERGENCE_ERROR
        case (-6)
            ! Pure relative error control requested but failed
            error stop DIFFEQ_ERROR_TEST_FAILURE
        end select

        ! Store the results
        if (nx == 2 .and. abs(t) > abs(xmax)) then
            call DVINDY(xmax, 0, rwork(21), neqn, y, istate)
            if (istate /= 0) error stop DIFFEQ_INVALID_OPERATION_ERROR
            t = xmax
        end if
        call this%append_to_buffer(t, y)

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
end subroutine

! ------------------------------------------------------------------------------
pure function vode_order_inquiry(this) result(rst)
    !! Gets the highest order of this integrator.  This integrator is a variable
    !! order integrator; therefore, the exact order is problem dependent.
    class(multistep_integrator), intent(in) :: this
        !! The multistep integrator.
    integer(int32) :: rst
        !! The highest order of this integrator.

    ! Process
    if (this%get_method() == DIFFEQ_ADAMS_METHOD) then
        rst = 12
    else
        rst = 5
    end if
end function

! ------------------------------------------------------------------------------
subroutine vode_eqn(neqn, t, y, ydot, rpar, ipar)
    !! The routine containing the differential equations solved by VODE.
    !!
    !! For a BDF problem with a supplied mass matrix, DVODE receives the
    !! equivalent explicit system \(y' = M^{-1}f(t,y)\).  The transformation
    !! is applied here so the public callback remains \(M y' = f\).
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
    if (associated(args%mass)) then
        call transform_mass_rhs(args%mass, neqn, t, y, ydot, args)
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine vode_jacobian(neqn, t, y, ml, mu, pd, nrowpd, rpar, ipar)
    !! The routine containing the Jacobian calculation routine.
    !!
    !! If a mass matrix is present, the supplied Jacobian is transformed to
    !! the Jacobian of the equivalent DVODE system, \(M^{-1}J\).  This matches
    !! the constant-mass formulation used by the BDF adapter.
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
    pd(1,1) = pd(1,1) + 0.0d0 * real(ml + mu, real64)
    if (associated(args%mass)) then
        call transform_mass_jacobian(args%mass, neqn, t, y, pd, args)
    end if
end subroutine

subroutine transform_mass_rhs(mass_callback, neqn, t, y, rhs, args)
    !! Transform a mass-matrix right-hand side into \(M^{-1}f\).
    procedure(ode_mass_matrix) :: mass_callback
    integer(int32), intent(in) :: neqn
    real(real64), intent(in) :: t, y(neqn)
    real(real64), intent(inout) :: rhs(neqn)
    class(vode_argument_container), intent(in) :: args
    real(real64) :: mass(neqn,neqn), solution(neqn)
    call mass_callback(t, y, mass, args%args)
    call solve_dense_system(mass, rhs, solution)
    rhs = solution
end subroutine

subroutine transform_mass_jacobian(mass_callback, neqn, t, y, jac, args)
    !! Transform Jacobian columns into the equivalent \(M^{-1}J\) system.
    procedure(ode_mass_matrix) :: mass_callback
    integer(int32), intent(in) :: neqn
    real(real64), intent(in) :: t, y(neqn)
    real(real64), intent(inout) :: jac(neqn,neqn)
    class(vode_argument_container), intent(in) :: args
    real(real64) :: mass(neqn,neqn), column(neqn)
    integer(int32) :: i
    call mass_callback(t, y, mass, args%args)
    do i = 1, neqn
        column = jac(:,i)
        call solve_dense_system(mass, column, jac(:,i))
    end do
end subroutine

subroutine solve_dense_system(matrix, rhs, solution)
    !! Solve a dense linear system with partial pivoting.
    real(real64), intent(in) :: matrix(:,:), rhs(:)
    real(real64), intent(out) :: solution(:)
    real(real64) :: a(size(rhs),size(rhs)), b(size(rhs)), factor, pivot_value
    real(real64) :: row(size(rhs))
    integer(int32) :: i, k, pivot, n
    n = size(rhs); a = matrix; b = rhs
    do k = 1, n - 1
        pivot = k; pivot_value = abs(a(k,k))
        do i = k + 1, n
            if (abs(a(i,k)) > pivot_value) then
                pivot = i; pivot_value = abs(a(i,k))
            end if
        end do
        if (pivot_value <= epsilon(1.0d0)) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
        if (pivot /= k) then
            row = a(k,:); a(k,:) = a(pivot,:); a(pivot,:) = row
            factor = b(k); b(k) = b(pivot); b(pivot) = factor
        end if
        do i = k + 1, n
            factor = a(i,k) / a(k,k)
            a(i,k:n) = a(i,k:n) - factor * a(k,k:n)
            b(i) = b(i) - factor * b(k)
        end do
    end do
    if (abs(a(n,n)) <= epsilon(1.0d0)) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    solution(n) = b(n) / a(n,n)
    do i = n - 1, 1, -1
        solution(i) = (b(i) - sum(a(i,i+1:n) * solution(i+1:n))) / a(i,i)
    end do
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

    rst = DIFFEQ_ADAMS_METHOD
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

    rst = DIFFEQ_BDF_METHOD
end function

! ------------------------------------------------------------------------------
end module