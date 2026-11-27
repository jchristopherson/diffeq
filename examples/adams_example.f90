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

program adams_example
    !! Integrate a non-stiff damped oscillator with the Adams PECE solver.
    !!
    !! The solver selects its step size and Adams order independently of the
    !! requested output times. The known solution provides a reference for
    !! both displacement and velocity.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use diffeq
    use fplot_core
    implicit none

    integer, parameter :: noutput = 401
    real(real64), parameter :: omega = 2.0d0, damping_ratio = 0.08d0
    real(real64), parameter :: final_time = 10.0d0
    real(real64), parameter :: damped_frequency = omega * sqrt(1.0d0 - damping_ratio**2)
    real(real64) :: times(noutput), exact_displacement(noutput), exact_velocity(noutput)
    real(real64) :: decay(noutput), max_error
    real(real64), allocatable :: solution(:,:)
    type(adams) :: integrator
    type(ode_container) :: model
    type(plot_2d) :: plt
    type(plot_data_2d) :: displacement, displacement_exact, velocity, velocity_exact
    class(legend), pointer :: key
    integer :: i, max_order_seen, steps

    ! Initialization
    do i = 1, noutput
        times(i) = final_time * real(i - 1, real64) / real(noutput - 1, real64)
    end do
    model%fcn => oscillator_rhs
    model%status_update => record_step
    max_order_seen = 0
    steps = 0
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)

    ! Compute the solution
    call integrator%solve(model, times, [1.0d0, 0.0d0])
    solution = integrator%get_solution()
    if (any(shape(solution) /= [noutput, 3])) error stop "Unexpected Adams output shape"
    if (.not.all(ieee_is_finite(solution))) error stop "Adams solution contains nonfinite values"

    ! Output
    decay = exp(-damping_ratio * omega * times)
    exact_displacement = decay * (cos(damped_frequency * times) + &
        damping_ratio * omega / damped_frequency * sin(damped_frequency * times))
    exact_velocity = -decay * omega**2 / damped_frequency * sin(damped_frequency * times)
    max_error = max(maxval(abs(solution(:,2) - exact_displacement)), &
        maxval(abs(solution(:,3) - exact_velocity)))
    print "(A,I0)", "Successful Adams steps: ", steps
    print "(A,I0)", "Highest Adams order used: ", max_order_seen
    print "(A,ES12.4)", "Maximum error against exact solution: ", max_error
    if (steps == 0 .or. max_order_seen == 0) error stop "No Adams steps were recorded"
    if (max_error > 1.0d-4) error stop "Adams solution failed the accuracy check"

    ! Plotting
    call plt%initialize()
    call plt%set_x_axis_title("t\n")
    call plt%set_y_axis_title("y(t) and dy/dt")
    key => plt%get_legend()
    call key%set_is_visible(.true.)
    call key%set_draw_border(.false.)
    call key%set_draw_inside_axes(.false.)
    call key%set_vertical_position(LEGEND_BOTTOM)
    call key%set_horizontal_position(LEGEND_CENTER)
    call key%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call displacement%define_data(solution(:,1), solution(:,2))
    call displacement%set_name("Displacement (Adams)")
    call displacement%set_line_width(2.0)
    call plt%push(displacement)
    call displacement_exact%define_data(times, exact_displacement)
    call displacement_exact%set_name("Displacement (exact)")
    call displacement_exact%set_line_style(LINE_DASHED)
    call plt%push(displacement_exact)
    call velocity%define_data(solution(:,1), solution(:,3))
    call velocity%set_name("Velocity (Adams)")
    call velocity%set_line_width(2.0)
    call plt%push(velocity)
    call velocity_exact%define_data(times, exact_velocity)
    call velocity_exact%set_name("Velocity (exact)")
    call velocity_exact%set_line_style(LINE_DASHED)
    call plt%push(velocity_exact)
    call plt%draw()

contains
    subroutine oscillator_rhs(t, y, dydt, args)
        !! Convert y'' + 2*zeta*omega*y' + omega**2*y = 0 to a first-order system.
        real(real64), intent(in) :: t, y(:)
        real(real64), intent(out) :: dydt(:)
        class(*), intent(in), optional :: args

        dydt(1) = y(2)
        dydt(2) = -omega**2 * y(1) - 2.0d0 * damping_ratio * omega * y(2)
    end subroutine

    subroutine record_step(info)
        !! Record the accepted-step count and largest order selected by Adams.
        type(ode_integration_status), intent(in) :: info

        steps = steps + 1
        max_order_seen = max(max_order_seen, info%order)
    end subroutine
end program adams_example