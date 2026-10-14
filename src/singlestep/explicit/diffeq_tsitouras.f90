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

module diffeq_tsitouras
    use iso_fortran_env
    use diffeq_base
    use diffeq_singlestep
    use diffeq_errors
    implicit none
    private
    public :: tsitouras_54

    type, extends(single_step_integrator) :: tsitouras_54
        !! The Tsitouras fifth-order Runge--Kutta method with an embedded
        !! fourth-order estimator.
        !!
        !! The seven-stage FSAL tableau is designed for efficient adaptive
        !! integration of non-stiff systems.
    contains
        procedure, public :: pre_step_action => tsit_pre_step
            !! Performs any pre-step actions.
        procedure, public :: get_order => tsit_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => tsit_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => tsit_get_stage_count
            !! Gets the stage count for this integrator.
        procedure, public :: attempt_step => tsit_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => tsit_post_step
            !! Performs any post-step actions.
        procedure, public :: interpolate => tsit_interpolate
            !! Performs the interpolation.
    end type

contains
! ------------------------------------------------------------------------------
subroutine tsit_pre_step(this, prevs, sys, h, x, y, f, args)
    !! Placeholder routine for any pre-step actions.
    class(tsitouras_54), intent(inout) :: this
        !! The tsitouras_54 object.
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

    ! Process
    return
end subroutine

! ------------------------------------------------------------------------------
pure function tsit_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(tsitouras_54), intent(in) :: this
        !! The tsitouras_54 object.
    integer :: rst
        !! The order.
    rst = 5
end function

! ------------------------------------------------------------------------------
pure function tsit_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(tsitouras_54), intent(in) :: this
        !! The tsitouras_54 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function tsit_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(tsitouras_54), intent(in) :: this
        !! The tsitouras_54 object.
    integer :: rst
        !! The stage count.
    rst = 7
end function

! ------------------------------------------------------------------------------
subroutine tsit_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    use diffeq_tsit45_constants
    !! Attempts an integration step for this integrator.
    !!
    !! The final stage is evaluated at the endpoint and is reused as the first
    !! stage of the next accepted step.  The fifth-order solution is paired
    !! with an embedded fourth-order solution; their difference supplies the
    !! local error estimate.
    class(tsitouras_54), intent(inout) :: this
        !! The tsitouras_54 object.
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
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Parameters
    integer, parameter :: nstages = 7

    ! Local Variables
    integer :: i, j
    real(real64) :: stage(size(y))

    ! Process
    ! k1 = f as the derivatives were computed from the previous step
    k(:,1) = f

    ! Compute the remaining stages
    do i = 2, nstages
        stage = y
        do j = 1, i - 1
            stage = stage + h * tsit_a(i,j) * k(:,j)
        end do
        call sys%fcn(x + tsit_c(i) * h, stage, k(:,i), args)
        call this%log_function_evaluation()
    end do

    ! Form the solution estimate and the error estimate
    yn = y
    yerr = 0.0d0
    do i = 1, nstages
        yn = yn + h * tsit_b(i) * k(:,i)
        yerr = yerr + h * tsit_e(i) * k(:,i)
    end do

    ! The FSAL stage supplies the derivatives at the end of the step
    call sys%fcn(x + h, yn, fn, args)
    call this%log_function_evaluation()
end subroutine

! ------------------------------------------------------------------------------
subroutine tsit_post_step(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    !! Placeholder routine for any post-step actions.
    !!
    !! The interpolation for this integrator requires only the solution and
    !! derivative values at each end of the step, so no additional storage is
    !! required here.
    class(tsitouras_54), intent(inout) :: this
        !! The tsitouras_54 object.
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
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Process
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine tsit_interpolate(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    !!
    !! A cubic Hermite polynomial is constructed from the solution and
    !! derivative values at each end of the step.
    class(tsitouras_54), intent(in) :: this
        !! The tsitouras_54 object.
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

    ! Local Variables
    real(real64) :: h, s

    ! Initialization
    h = xn1 - xn
    s = (x - xn) / h

    ! Process
    y = (2.0d0 * s**3 - 3.0d0 * s**2 + 1.0d0) * yn + &
        (s**3 - 2.0d0 * s**2 + s) * h * fn + &
        (-2.0d0 * s**3 + 3.0d0 * s**2) * yn1 + &
        (s**3 - s**2) * h * fn1
end subroutine

! ------------------------------------------------------------------------------
end module