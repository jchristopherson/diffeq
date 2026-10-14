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

program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: mass = 2.0d0

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable, dimension(:) :: theta
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: pltxy, plttheta
    class(legend), pointer :: lgnd

    ! Pass in the model parameters
    args%length = length
    args%mass = mass

    ! Define the model
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)

    ! Compute the solution
    call integrator%solve( &
        mdl, &                                  ! model to solve
        [0.0d0, 1.0d1], &                       ! time bounds
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! initial conditions - must satisfy algebraic constraints
        args = args &                           ! arguments to pass to the model
    )
    sol = integrator%get_solution()

    ! Compute the angular representation of the solution
    theta = atan2(sol(:,4), sol(:,2)) * 1.8d2 / pi

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
    call pltxy%push(sol(:,1), sol(:,4), name = "y(t)")

    call plt%set(1, 1, pltxy)

    ! Plot the angular results
    call plttheta%initialize()
    call plttheta%set_x_axis_title("t")
    call plttheta%set_y_axis_title("{/Symbol q} [deg]")
    call plttheta%push(sol(:,1), theta)
    
    call plt%set(1, 2, plttheta)

    call plt%draw()
end program