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

module diffeq_c_integrators
    use iso_fortran_env
    use iso_c_binding
    use diffeq_c_containers
    use diffeq_c_errors
    use diffeq
    implicit none
    private
    public :: c_ode
    public :: c_ode_jacobian
    public :: c_ode_mass_matrix
    public :: c_ode_settings
    public :: c_default_ode_settings

    integer, parameter :: C_RUNGE_KUTTA_23 = 1
    integer, parameter :: C_RUNGE_KUTTA_45 = 2
    integer, parameter :: C_RUNGE_KUTTA_853 = 3
    integer, parameter :: C_TSITOURAS = 4
    integer, parameter :: C_ROSENBROCK = 5
    integer, parameter :: C_KCARP4 = 6
    integer, parameter :: C_KCARP5 = 7
    integer, parameter :: C_BDF = 8
    integer, parameter :: C_ADAMS = 9

    interface
        subroutine c_ode(n, x, y, dydx)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: dydx(n)
        end subroutine

        subroutine c_ode_jacobian(m, n, x, y, jac, ldj)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: m
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: jac(ldj,n)
            integer(c_int), intent(in), value :: ldj
        end subroutine

        subroutine c_ode_mass_matrix(n, x, y, m, ldm)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: m(ldm,n)
            integer(c_int), intent(in), value :: ldm
        end subroutine
    end interface

    type, bind(C) :: c_ode_settings
        !! A container for ODE settings.
        logical(c_bool) :: mass_matrix_constant
            !! True if the mass matrix is constant regardless of state, or
            !! false if the mass matrix is state-dependent.
        real(c_double) :: minimum_step_size
            !! The minimum allowed integrator step size.
        real(c_double) :: maximum_step_size
            !! The maximum allowed integrator step size.
        integer(c_int) :: step_limit
            !! The limit on the maximum total number of integration steps.
        logical(c_bool) :: allow_overshoot
            !! True if the solver is allowed to overshoot the final value in
            !! the integration range and interpolate back; else, false if the
            !! solver must terminate on the final value.
        real(c_double) :: beta
            !! PI step size controller exponent.  Set to zero to turn off PI
            !! step size control.
        real(c_double) :: absolute_tolerance
            !! The absolute tolerance applied to each equation.
        real(c_double) :: relative_tolerance
            !! The relative tolerance applied to each equation.
    end type

contains
! ------------------------------------------------------------------------------
    subroutine c_default_ode_settings(settings) &
        bind(C, name = "c_default_ode_settings")
        !! Fills the supplied [[c_ode_settings]] object with default settings.
        type(c_ode_settings), intent(inout) :: settings
            !! The [[c_ode_settings]] object.
        settings%absolute_tolerance = 1.0d-6
        settings%allow_overshoot = .true.
        settings%beta = 0.0d0
        settings%mass_matrix_constant = .true.
        settings%maximum_step_size = huge(1.0d0)
        settings%minimum_step_size = 1.0d1 * epsilon(1.0d0)
        settings%relative_tolerance = 1.0d-6
        settings%step_limit = 1000000
    end subroutine

! ------------------------------------------------------------------------------
    function c_integrate_ode(method, nx, nstates, fcn, x, initcond, &
        settings, jac, mass) result(rst)
        integer(c_int), intent(in) :: method
        integer(c_int), intent(in) :: nx
        integer(c_int), intent(in) :: nstates
        procedure(c_ode), intent(in), pointer :: fcn
        real(c_double), intent(in) :: x(nx)
        real(c_double), intent(in) :: initcond(nstates)
        type(c_ode_settings), intent(in) :: settings
        procedure(c_ode_jacobian), intent(in), pointer :: jac
        procedure(c_ode_mass_matrix), intent(in), pointer :: mass
        integer(c_int) :: rst

        class(ode_integrator), pointer :: integrator
        type(runge_kutta_23), target :: rk23
        type(runge_kutta_45), target :: rk45
        type(runge_kutta_853), target :: rk853
        type(tsitouras_54), target :: tsit
        type(rosenbrock), target :: rsbrk
        type(kennedy_carpenter_4), target :: kc4
        type(kennedy_carpenter_5), target :: kc5
        type(bdf), target :: bdf_i
        type(adams), target :: adams_i
        type(ode_container) :: model

        rst = C_NO_ERROR
        select case (method)
        case (C_RUNGE_KUTTA_23)
            integrator => rk23
        case (C_RUNGE_KUTTA_45)
            integrator => rk45
        case (C_RUNGE_KUTTA_853)
            integrator => rk853
        case (C_TSITOURAS)
            integrator => tsit
        case (C_ROSENBROCK)
            integrator => rsbrk
        case (C_KCARP4)
            integrator => kc4
        case (C_KCARP5)
            integrator => kc5
        case (C_BDF)
            integrator => bdf_i
        case (C_ADAMS)
            integrator => adams_i
        case default
            rst = C_INVALID_INPUT_ERROR
            return
        end select

    contains
        subroutine internal_ode(x_, y_, dydx_, args_)
            real(real64), intent(in) :: x_
            real(real64), intent(in), dimension(:) :: y_
            real(real64), intent(out), dimension(:) :: dydx_
            class(*), intent(inout), optional :: args_
            call fcn(size(y_), x_, y_, dydx_)
        end subroutine
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module