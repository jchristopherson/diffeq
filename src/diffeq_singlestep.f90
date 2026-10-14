! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module diffeq_singlestep
    use iso_fortran_env
    use diffeq_base
    use diffeq_errors
    implicit none
    private
    public :: attempt_single_step
    public :: get_single_step_logical_parameter
    public :: single_step_post_step_routine
    public :: single_step_pre_step_routine
    public :: single_step_interpolate
    public :: single_step_integer_inquiry
    public :: single_step_integrator
    
    ! ------------------------------------------------------------------------------
    type, abstract, extends(ode_integrator) :: single_step_integrator
        !! The most basic, single-step integrator object capable of integrating
        !! systems of ODE's.
    contains
        procedure(attempt_single_step), public, pass, deferred :: attempt_step
            !! Attempts an integration step for a single-step integrator.
        procedure(get_single_step_logical_parameter), public, pass, &
            deferred :: get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure(single_step_post_step_routine), public, pass, deferred :: &
            post_step_action
            !! Performs actions such as setting up interpolation after 
            !! completion of a successful integration step.
        procedure(single_step_interpolate), public, pass, deferred :: &
            interpolate
            !! Performs an interpolation to estimate the solution at the
            !! requested point.
        procedure(single_step_pre_step_routine), public, pass, deferred :: &
            pre_step_action
            !! Provides a routine for performing any actions, such as setting
            !! up Jacobian calculations.
        procedure(single_step_integer_inquiry), public, pass, deferred :: &
            get_stage_count
            !! Gets the number of stages used by the integrator.
        procedure, public :: solve => ssi_ode_solver
            !! Solves the supplied system of ODE's.
    end type

    interface
        subroutine attempt_single_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
            use iso_fortran_env
            import single_step_integrator
            import ode_container
            !! Attempts an integration step for a single-step integrator.
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
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
                !! An N-by-NSTAGES matrix containing the derivatives at each
                !! stage.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        pure function get_single_step_logical_parameter(this) result(rst)
            !! Returns a logical parameter from a single_step_integrator object.
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
            logical :: rst
                !! The parameter.
        end function

        subroutine single_step_post_step_routine(this, sys, dense, x, xn, y, &
            yn, f, fn, k, args)
            !! Provides a routine for performing any actions, such as setting
            !! up interpolation, after successful completion of a step.
            use iso_fortran_env
            import single_step_integrator
            import ode_container
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
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
                !! An N-by-NSTAGES matrix containing the derivatives at each
                !! stage.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        subroutine single_step_interpolate(this, x, xn, yn, fn, xn1, yn1, &
            fn1, y)
            !! Provides a routine for interpolation.
            use iso_fortran_env
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
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
        end subroutine

        subroutine single_step_pre_step_routine(this, prevs, sys, h, x, y, f, &
            args)
            !! Provides a routine for performing any actions, such as setting
            !! up Jacobian calculations.
            use iso_fortran_env
            import single_step_integrator
            import ode_container
            class(single_step_integrator), intent(inout) :: this
                !! The single_step_integrator object.
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
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the differential equation subroutine.
        end subroutine

        pure function single_step_integer_inquiry(this) result(rst)
            !! Gets an integer from the integrator.
            use iso_fortran_env
            import single_step_integrator
            class(single_step_integrator), intent(in) :: this
                !! The single_step_integrator object.
            integer :: rst
                !! The integer value.
        end function
    end interface

contains
! ******************************************************************************
! SINGLE-STEP INTEGRATOR
! ------------------------------------------------------------------------------
subroutine ssi_ode_solver(this, sys, x, iv, args)
    !! Solves the supplied system of ODE's.
    class(single_step_integrator), intent(inout) :: this
        !! The single_step_integrator object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in), dimension(:) :: x
        !! An array of independent variable values at which to return the 
        !! the solution to the ODE's.
    real(real64), intent(in), dimension(:) :: iv
        !! An array containing the initial values for each ODE.
    class(*), intent(inout), optional, target :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    logical :: dense, success, forward
    integer :: i, j, n, neqn, nsteps, nstages
    real(real64) :: h, xo, xn, xmax, ei, eold
    real(real64), allocatable, dimension(:) :: f, y, yn, fn, yerr, yi
    real(real64), allocatable, dimension(:,:) :: k
    
    ! Initialization
    n = size(x)
    success = .true.

    ! Input Checking
    if (n < 2) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (.not.sys%get_is_ode_defined()) error stop DIFFEQ_MISSING_ARGUMENT_ERROR

    ! Additional Initialization
    neqn = size(iv)
    xo = x(1)
    xmax = x(n)
    forward = xmax >= xo
    j = 2
    eold = 1.0d-4
    dense = (n > 2)
    nsteps = this%get_step_limit()
    nstages = this%get_stage_count()
    
    ! Memory Allocations
    allocate( &
        f(neqn), &
        y(neqn), &
        yn(neqn), &
        fn(neqn), &
        yerr(neqn), &
        k(neqn, nstages), &
        yi(neqn) &
    )

    ! Estimate an initial step size
    !
    ! Outputs:
    ! - f: Value of the derivatives at xo
    ! - h: Initial step size estimate
    call this%estimate_inital_step_size(sys, xo, xmax, iv, f, h, args)
    
    ! Store the initial conditions
    call this%append_to_buffer(x(1), iv)
    y = iv

    ! Cycle until integration is complete
    do i = 1, nsteps
        ! Perform any pre-step actions
        call this%pre_step_action(success, sys, h, xo, y, f, args = args)
        
        ! Attempt a step
        call this%attempt_step(sys, h, xo, y, f, yn, fn, yerr, k, args)
        xn = xo + h

        ! Compute a normalized error value.  A value < 1 indicates success
        ei = this%compute_error_norm(y, yn, yerr)

        ! Determine the next step size
        h = this%estimate_next_step_size(ei, eold, h, xo)

        ! Reject the step?
        success = ei <= 1.0d0
        if (.not.success) cycle   ! We failed.  Try again with a new step size

        ! If we're here, the step has been successful.  Take any post-step
        ! action such as setting up interpolation routines, etc.
        call this%post_step_action(sys, dense, xo, xn, y, yn, f, fn, k, args)

        ! Do we need to interpolate for dense output, or can we just store
        ! values and move on
        if (dense) then
            ! Perform the interpolation as needed until a new step is required
            interp : do while (j <= n .and. &
                ((forward .and. x(j) <= xn) .or. &
                (.not.forward .and. x(j) >= xn)))
                call this%interpolate(x(j), xo, y, f, xn, yn, fn, yi)
                call this%append_to_buffer(x(j), yi)
                j = j + 1
                if (j > n) exit interp
            end do interp
        else
            ! Store the values and move on
            call this%append_to_buffer(xn, yn)
        end if

        ! Are we done?
        if ((xmax - x(1)) * (xn - xmax) >= 0.0d0) then
            ! Deal with the case where output is only returned at the
            ! integration points and the solver oversteps the endpoint.
            if ((xmax - x(1)) * (xn - xmax) > 0.0d0) then
                ! Interpolate to get the solution at xmax
                call this%post_step_action(sys, .true., xo, xn, y, yn, f, fn, k)
                if (j <= n) then
                    call this%interpolate(xmax, xo, y, f, xn, yn, fn, yi)
                    call this%append_to_buffer(xmax, yi)
                end if
            end if
            
            ! We're done
            go to 100
        end if

        ! Do we need to limit the step size to not overshoot the terminal value?
        if (this%get_allow_overshoot() .and. &
            (xmax - x(1)) * (xn + h - xmax) > 0.0d0) then
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
            call sys%fcn(xn, yn, f, args) ! write to f, not fn
        end if
    end do

    ! If we're here, the solver has run out of allowable steps
    error stop DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR

    ! End
100 continue
    return
end subroutine

! ------------------------------------------------------------------------------
end module