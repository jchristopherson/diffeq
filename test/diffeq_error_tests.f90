! This file is part of diffeq.
! 
! diffeq is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! diffeq is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with diffeq. If not, see <https://www.gnu.org/licenses/>.

program diffeq_error_tests
    use diffeq
    implicit none

    character(len=32) :: test_name
    type(runge_kutta_45) :: integrator
    type(ode_container) :: mdl

    call get_command_argument(1, test_name)
    select case (trim(test_name))
    case ("missing_ode")
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("short_grid")
        mdl%fcn => constant_ode
        call integrator%solve(mdl, [0.0d0], [0.0d0])
    case ("step_limit")
        mdl%fcn => constant_ode
        call integrator%set_step_limit(0)
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("tolerance_size")
        mdl%fcn => constant_ode
        call integrator%set_absolute_tolerance([1.0d-6, 1.0d-6])
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case default
        error stop 2
    end select

contains
    subroutine constant_ode(x, y, dydx, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: dydx(:)
        class(*), intent(inout), optional :: args
        dydx = 1.0d0
    end subroutine
end program