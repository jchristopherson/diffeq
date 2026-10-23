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

module diffeq_rk45
    use iso_fortran_env
    use diffeq_base
    use diffeq_singlestep
    use diffeq_errors
    implicit none
    private
    public :: runge_kutta_45

    type, extends(single_step_integrator) :: runge_kutta_45
        !! The Dormand--Prince RK45 integrator.
        !!
        !! A fifth-order solution is paired with an embedded fourth-order
        !! solution.  The estimated local defect is
        !! \f$ e_{n+1} = y_{n+1}^{[5]} - y_{n+1}^{[4]} \f$, and the method uses
        !! the fifth-order result after an accepted step.  The method has
        !! seven stages and FSAL (first-same-as-last) reuse, so the final
        !! stage derivative becomes the first derivative of the next step.
        !! It is a general-purpose choice for moderate-to-high accuracy on
        !! non-stiff problems.
        real(real64), private, allocatable, dimension(:) :: rc1
        real(real64), private, allocatable, dimension(:) :: rc2
        real(real64), private, allocatable, dimension(:) :: rc3
        real(real64), private, allocatable, dimension(:) :: rc4
        real(real64), private, allocatable, dimension(:) :: rc5
    contains
        procedure, public :: pre_step_action => rk45_pre_step
            !! Performs any pre-step actions.
        procedure, public :: get_order => rk45_get_order
            !! Gets the order of the integrator.
        procedure, public :: get_is_fsal => rk45_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => rk45_get_stage_count
            !! Gets the stage count for this integrator.
        procedure, private :: initialize_interp => rk45_init_interp
            !! Allocates storage for the interpolation process.
        procedure, public :: attempt_step => rk45_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => rk45_set_up_interp
            !! Sets up the interpolation process as the post-step action.
        procedure, public :: interpolate => rk45_interp
            !! Performs the interpolation.
    end type

contains
! ------------------------------------------------------------------------------
subroutine rk45_pre_step(this, prevs, sys, h, x, y, f, args)
    !! Placeholder routine for any pre-step actions.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
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
pure function rk45_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    integer :: rst
        !! The order.
    rst = 5
end function

! ------------------------------------------------------------------------------
pure function rk45_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure function rk45_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
    integer :: rst
        !! The stage count.
    rst = 7
end function

! ------------------------------------------------------------------------------
subroutine rk45_init_interp(this, neqn)
    !! Allocates storage for the interpolation process.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
    integer, intent(in) :: neqn
        !! The number of equations being integrated.

    ! Process
    if (allocated(this%rc1)) then
        if (size(this%rc1) == neqn) then
            ! All is good
            return
        else
            deallocate(this%rc1)
            deallocate(this%rc2)
            deallocate(this%rc3)
            deallocate(this%rc4)
            deallocate(this%rc5)
        end if
    end if

    allocate( &
        this%rc1(neqn), &
        this%rc2(neqn), &
        this%rc3(neqn), &
        this%rc4(neqn), &
        this%rc5(neqn) &
    )
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    use diffeq_dprk45_constants
    !! Attempts an integration step for this integrator.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
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

    ! Local Variables
    integer :: n

    ! Initialization
    n = size(y)

    ! Process
    ! k1 = f as the derivatives were computed from the previous step
    k(:,1) = f

    yn = y + h * a21 * f
    call sys%compute_rhs(x + h * c2, yn, k(:,2), args)

    yn = y + h * (a31 * f + a32 * k(:,2))
    call sys%compute_rhs(x + h * c3, yn, k(:,3), args)

    yn = y + h * (a41 * f + a42 * k(:,2) + a43 * k(:,3))
    call sys%compute_rhs(x + h * c4, yn, k(:,4), args)

    yn = y + h * (a51 * f + a52 * k(:,2) + a53 * k(:,3) + a54 * k(:,4))
    call sys%compute_rhs(x + h * c5, yn, k(:,5), args)

    yn = y + h * (a61 * f + a62 * k(:,2) + a63 * k(:,3) + a64 * k(:,4) + &
        a65 * k(:,5))
    call sys%compute_rhs(x + h * c6, yn, k(:,6), args)

    yn = y + h * (a71 * f + a73 * k(:,3) + a74 * k(:,4) + a75 * k(:,5) + a76 * k(:,6))
    call sys%compute_rhs(x + h * c7, yn, fn, args)

    call this%log_function_evaluation(6)

    ! Compute the error estimate
    yerr = h * (e1 * f + e2 * k(:,2) + e3 * k(:,3) + e4 * k(:,4) + &
        e5 * k(:,5) + e6 * k(:,6) + e7 * fn)
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_set_up_interp(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    use diffeq_dprk45_constants
    !! Sets up the interpolation process.
    class(runge_kutta_45), intent(inout) :: this
        !! The runge_kutta_45 object.
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

    ! Local Variables
    integer :: i, n
    real(real64) :: h, ydiff, bspl

    ! Quick Return
    if (.not.dense) return

    ! Initialization
    n = size(y)
    call this%initialize_interp(n)
    h = xn - x

    ! Process
    do i = 1, n
        this%rc1(i) = y(i)
        ydiff = yn(i) - y(i)
        this%rc2(i) = ydiff
        bspl = h * f(i) - ydiff
        this%rc3(i) = bspl
        this%rc4(i) = ydiff - h * fn(i) - bspl
        this%rc5 = h * (d1 * f(i) + d3 * k(i,3) + d4 * k(i,4) + &
            d5 * k(i,5) + d6 * k(i,6) + d7 * fn(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine rk45_interp(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    class(runge_kutta_45), intent(in) :: this
        !! The runge_kutta_45 object.
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
    real(real64) :: h, s, s1

    ! Process
    h = xn1 - xn
    s = (x - xn) / h
    s1 = 1.0d0 - s
    y = this%rc1 + s * (this%rc2 + s1 * (this%rc3 + &
        s * (this%rc4 + s1 * this%rc5)))
end subroutine

! ------------------------------------------------------------------------------
end module