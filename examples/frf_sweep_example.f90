module ode_module
    use iso_fortran_env
    use diffeq_harmonics
    use ieee_arithmetic
    implicit none

    ! Model Constants
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 0.04d0
    real(real64), parameter :: delta = 0.1d0
    real(real64), parameter :: gamma = 1.0d0

    ! ODE_CONTAINER object
    type, extends(harmonic_ode_container) :: duffing_ode
    contains
        ! Overloads the ODE subroutine providing the ability to obtain
        ! frequency information from the solver via the excitation_frequency
        ! property.
        procedure, public :: ode => duffing_model
    end type

contains
    subroutine duffing_model(this, x, y, dydx)
        ! Arguments
        class(duffing_ode), intent(in) :: this
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: dydx

        ! Equations
        dydx(1) = y(2)
        dydx(2) = gamma * cos(this%excitation_frequency * x) - &
            delta * y(2) - alpha * y(1) - beta * y(1)**3
    end subroutine

    ! Analytical Solution - Leg 1
    pure elemental function leg1(z) result(rst)
        ! Arguments
        real(real64), intent(in) :: z
        real(real64) :: rst

        ! Local Variables
        real(real64) :: arg, s, nan

        ! Process
        nan = ieee_value(nan, ieee_quiet_nan)
        arg = 4.0d0 * gamma**2 - 3.0d0 * beta * delta**2 * z**4 + &
            (delta**2 - 4.0d0 * alpha) * delta**2 * z**2
        if (arg < 0.0d0) then
            s = -1.0d0
        else
            s = (2.0d0 * sqrt(arg) + z * (3.0d0 * beta * z**2 - &
                2.0d0 * delta**2 + 4.0d0 * alpha)) / (4.0d0 * z)
        end if
        if (s < 0.0d0) then
            rst = nan
        else
            rst = sqrt(s)
        end if
    end function

    ! Analytical Solution - Leg 2
    pure elemental function leg2(z) result(rst)
        ! Arguments
        real(real64), intent(in) :: z
        real(real64) :: rst

        ! Local Variables
        real(real64) :: arg, s, nan

        ! Process
        nan = ieee_value(nan, ieee_quiet_nan)
        arg = 4.0d0 * gamma**2 - 3.0d0 * beta * delta**2 * z**4 + &
            (delta**2 - 4.0d0 * alpha) * delta**2 * z**2
        if (arg < 0.0d0) then
            s = -1.0d0
        else
            s = (-2.0d0 * sqrt(arg) + z * (3.0d0 * beta * z**2 - &
                2.0d0 * delta**2 + 4.0d0 * alpha)) / (4.0d0 * z)
        end if
        if (s < 0.0d0) then
            rst = nan
        else
            rst = sqrt(s)
        end if
    end function
end module

program example
    use iso_fortran_env
    use ode_module
    use diffeq_harmonics
    use fplot_core
    implicit none
    
    ! Parameters
    real(real64), parameter :: f1 = 0.5d0
    real(real64), parameter :: f2 = 2.0d0
    real(real64), parameter :: df = 0.05d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: deg = 1.8d2 / pi

    ! Local Variables
    type(duffing_ode) :: sys
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:) :: fup, fdown, w1, w2, z
    complex(real64), allocatable, dimension(:,:) :: rup, rdown

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    type(plot_data_2d) :: pd1, pd2, pd3
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the frequency vectors
    n = floor((f2 - f1) / df) + 1
    allocate(fup(n), fdown(n))
    fup = (/ (df * i + f1, i = 0, n - 1) /)
    fdown = (/ (df * i + f1, i = n - 1, 0, -1) /)

    ! Perform the ascending sweep
    rup = frequency_response(sys, fup, [0.0d0, 0.0d0])

    ! Perform the descending sweep
    rdown = frequency_response(sys, fdown, [0.0d0, 0.0d0])

    ! Compute the analytical solution
    z = linspace(0.5d0, 1.0d1, 1000)
    w1 = leg1(z)
    w2 = leg2(z)

    ! Plot the results
    call plt%initialize(2, 1)
    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call xAxis%set_title("w")
    call yAxis%set_title("X(w)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_LEFT)

    call pd1%define_data(fup, abs(rup(:,1)))
    call pd1%set_name("Ascending")
    call pd1%set_line_width(2.0)
    call pd1%set_draw_line(.false.)
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_style(MARKER_EMPTY_CIRCLE)
    call plt1%push(pd1)

    call pd2%define_data(fdown, abs(rdown(:,1)))
    call pd2%set_name("Descending")
    call pd2%set_line_width(2.0)
    call pd2%set_draw_line(.false.)
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_EMPTY_TRIANGLE)
    call plt1%push(pd2)

    call pd3%define_data(w1, z)
    call pd3%set_name("Analytical")
    call pd3%set_line_width(2.0)
    call pd3%set_line_color(CLR_BLACK)
    call plt1%push(pd3)

    call pd3%define_data(w2, z)
    call pd3%set_name("")
    call plt1%push(pd3)

    call plt2%initialize()
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("w")
    call yAxis%set_title("{/Symbol f} [deg]")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 2.0d0)

    call pd1%define_data(fup, deg * atan2(aimag(rup(:,1)), real(rup(:,1))))
    call plt2%push(pd1)

    call pd2%define_data(fdown, deg * atan2(aimag(rdown(:,1)), real(rdown(:,1))))
    call plt2%push(pd2)

    call plt%set(1, 1, plt1)
    call plt%set(2, 1, plt2)
    call plt%draw()
end program