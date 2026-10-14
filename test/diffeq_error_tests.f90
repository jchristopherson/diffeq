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
    type(bdf) :: bdf_integrator
    type(adams) :: adams_integrator
    type(ode_container) :: mdl
    type(ggl_container) :: ggl_mdl
    real(real64), dimension(8) :: ggl_state, ggl_rhs
    real(real64), dimension(8,8) :: ggl_mass
    logical :: consistent

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
    case ("bdf_step_exhaustion")
        mdl%fcn => constant_ode
        call bdf_integrator%set_step_limit(1)
        call bdf_integrator%set_maximum_step_size(1.0d-3)
        call bdf_integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("adams_step_exhaustion")
        mdl%fcn => constant_ode
        call adams_integrator%set_step_limit(1)
        call adams_integrator%set_maximum_step_size(1.0d-3)
        call adams_integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("ggl_dimensions")
        call ggl_mdl%initialize(2, 3)
    case ("ggl_gains")
        call ggl_mdl%initialize(2, 1, alpha = -1.0d0)
    case ("ggl_state_size", "ggl_rhs_size", "ggl_matrix_size", "ggl_missing_callback", &
        "ggl_tolerance", "ggl_singular_mass", "ggl_singular_constraints")
        call ggl_mdl%initialize(2, 1)
        ggl_mdl%fcn => constant_ode
        ggl_mdl%constraints => linear_constraints
        ggl_mdl%constraint_jacobian => linear_gradient
        ggl_mdl%constraint_curvature => linear_curvature
        ggl_state = 0.0d0
        select case (trim(test_name))
        case ("ggl_state_size")
            call ggl_mdl%compute_rhs(0.0d0, ggl_state(1:5), ggl_rhs(1:6))
        case ("ggl_rhs_size")
            call ggl_mdl%compute_rhs(0.0d0, ggl_state(1:6), ggl_rhs(1:5))
        case ("ggl_matrix_size")
            call ggl_mdl%compute_mass_matrix(0.0d0, ggl_state(1:6), ggl_mass(1:5,1:5))
        case ("ggl_missing_callback")
            nullify(ggl_mdl%constraint_curvature)
            call ggl_mdl%compute_rhs(0.0d0, ggl_state(1:6), ggl_rhs(1:6))
        case ("ggl_tolerance")
            consistent = ggl_mdl%validate_initial_conditions(0.0d0, ggl_state(1:6), -1.0d0)
        case ("ggl_singular_mass")
            ggl_mdl%mass_matrix => singular_mass
            call ggl_mdl%compute_rhs(0.0d0, ggl_state(1:6), ggl_rhs(1:6))
        case ("ggl_singular_constraints")
            call ggl_mdl%initialize(2, 2)
            call ggl_mdl%initialize_multipliers(0.0d0, ggl_state)
        end select
    case default
        error stop 2
    end select

contains
    subroutine linear_constraints(x, q, g, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(inout), optional :: args
        g = q(1)
    end subroutine

    subroutine linear_gradient(x, q, gradient, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: gradient
        class(*), intent(inout), optional :: args
        gradient = 0.0d0
        gradient(:,1) = 1.0d0
    end subroutine

    subroutine linear_curvature(x, q, direction, v, curvature, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q, direction, v
        real(real64), intent(out), dimension(:) :: curvature
        class(*), intent(inout), optional :: args
        curvature = 0.0d0
    end subroutine

    subroutine singular_mass(x, q, mass, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: mass
        class(*), intent(inout), optional :: args
        mass = 0.0d0
    end subroutine

    subroutine constant_ode(x, y, dydx, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: dydx(:)
        class(*), intent(inout), optional :: args
        dydx = 1.0d0
    end subroutine
end program