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

program stabilized_pendulum
    use iso_fortran_env, only : real64
    use diffeq, only : ggl_container, bdf
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: body_mass = 2.0d0

    ! Local Variables
    type(ggl_container) :: model
    type(bdf) :: solver
    real(real64), dimension(6) :: initial_state
    real(real64), allocatable, dimension(:) :: theta
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: pltxy, plttheta
    class(legend), pointer :: lgnd

    ! Set up the model and constraints
    call model%initialize(2, 1, alpha = 2.0d0, beta = 3.0d0)
    model%fcn => force
    model%mass_matrix => physical_mass
    model%constraints => position_constraint
    model%constraint_jacobian => constraint_gradient
    model%constraint_curvature => constraint_curvature
    call model%set_is_mass_matrix_dependent(.false.)

    ! Define the initial conditions and validate against the constraints.  The
    ! initial conditions must satisfy the constraints.
    initial_state = [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    if (.not.model%validate_initial_conditions(0.0d0, initial_state, 1.0d-10)) &
        error stop "Inconsistent physical initial conditions"
    call model%initialize_multipliers(0.0d0, initial_state)

    ! Set up the solver and adjust tolerances.  Utilize the same tolerance for
    ! every equation.
    call solver%set_absolute_tolerance(1.0d-9)
    call solver%set_relative_tolerance(1.0d-9)

    ! Solve the model and retrieve the results
    call solver%solve(model, [0.0d0, 1.0d1], initial_state)
    sol = solver%get_solution()

    ! Compute the angular representation of the solution
    theta = atan2(sol(:,3), sol(:,2)) * 1.8d2 / pi

! ------------------------------------------------------------------------------
    ! Set up the plot
    call plt%initialize(1, 2, width = 1200)

    ! Plot the X-Y results
    call pltxy%initialize()
    lgnd => pltxy%get_legend()
    call pltxy%set_x_axis_title("t")
    call pltxy%set_y_axis_title("x(t), y(t)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pltxy%push(sol(:,1), sol(:,2), name = "x(t)")
    call pltxy%push(sol(:,1), sol(:,3), name = "y(t)")

    call plt%set(1, 1, pltxy)

    ! Plot the angular results
    call plttheta%initialize()
    call plttheta%set_x_axis_title("t")
    call plttheta%set_y_axis_title("{/Symbol q} [deg]")
    call plttheta%push(sol(:,1), theta)
    
    call plt%set(1, 2, plttheta)

    call plt%draw()
contains
    subroutine force(x, y, applied_force, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: applied_force
        class(*), intent(inout), optional :: args
        applied_force = [0.0d0, -body_mass * 9.81d0]
    end subroutine

    subroutine physical_mass(x, q, mass, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: mass
        class(*), intent(inout), optional :: args
        mass = 0.0d0
        mass(1,1) = body_mass
        mass(2,2) = body_mass
    end subroutine

    subroutine position_constraint(x, q, g, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(inout), optional :: args
        g(1) = dot_product(q, q) - length**2
    end subroutine

    subroutine constraint_gradient(x, q, gradient, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: gradient
        class(*), intent(inout), optional :: args
        gradient(1,:) = 2.0d0 * q
    end subroutine

    subroutine constraint_curvature(x, q, direction, v, curvature, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q, direction, v
        real(real64), intent(out), dimension(:) :: curvature
        class(*), intent(inout), optional :: args
        curvature(1) = 2.0d0 * dot_product(direction, v)
    end subroutine
end program
