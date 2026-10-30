program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    integer, parameter :: npts = 1000
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]
    ! Timing Variables
    integer(int64) :: tsAD, tfAD, tsFD, tfFD, ts, tf, rate
    real(real64) :: timeAD, timeFD, time, drate, x(npts)

    ! Variables
    type(ode_container) :: modelAD, modelFD, model
    type(bdf) :: integrator
    type(ode_results) :: rstAD, rstFD, rst

    ! Set up the timing
    call system_clock(count_rate = rate)
    drate = real(rate, real64)

    ! Set up the solution points for a dense solve
    x = linspace(0.0d0, 5.0d1, npts)

    ! Set up the AD version and solve
    modelAD%fcn_ad => vanderpol_ad
    call system_clock(tsAD)
    call integrator%solve(modelAD, x, ic)
    call system_clock(tfAD)
    timeAD = real(tfAD - tsAD, real64) / drate
    rstAD = integrator%get_ode_results()
    print "(A)", "AD RESULTS:"
    call print_statistics(rstAD, timeAD)

    ! Set up the FD version and solve
    modelFD%fcn => vanderpol
    call integrator%clear_buffer()
    call system_clock(tsFD)
    call integrator%solve(modelFD, x, ic)
    call system_clock(tfFD)
    timeFD = real(tfFD - tsFD, real64) / drate
    rstFD = integrator%get_ode_results()
    print "(A)", "FD RESULTS:"
    call print_statistics(rstFD, timeFD)

    ! Set up the analytical Jacobian version
    model%jacobian => vanderpol_analytical_jacobian
    model%fcn => vanderpol
    call integrator%clear_buffer()
    call system_clock(ts)
    call integrator%solve(model, x, ic)
    call system_clock(tf)
    time = real(tf - ts, real64) / drate
    rst = integrator%get_ode_results()
    print "(A)", "ANALYTICAL JACOBIAN RESULTS:"
    call print_statistics(rst, time)

    ! Plot the results
    call plot_results()
contains

subroutine vanderpol_ad(x_, y_, dydx_, args_)
    type(real_ad), intent(in) :: x_, y_(:)
    type(real_ad), intent(out) :: dydx_(:)
    class(*), intent(inout), optional :: args_

    ! Model Constants
    real(real64), parameter :: mu_ = 5.0d0

    ! Equations
    dydx_(1) = y_(2)
    dydx_(2) = mu_ * (1.0d0 - y_(1)**2) * y_(2) - y_(1)
end subroutine

subroutine print_statistics(rst_, t_)
    type(ode_results), intent(in) :: rst_
    real(real64), intent(in) :: t_
    print "(A, F0.2, A)", achar(9) // "Execution Time: ", 1.0d3 * t_, " ms"
    print "(A, I0)", achar(9) // "Successful Steps: ", rst_%statistics%successful_steps
    print "(A, I0)", achar(9) // "Failed Steps: ", rst_%statistics%failed_steps
    print "(A, I0)", achar(9) // "Jacobian Evaluations: ", rst_%statistics%jacobian_evaluations
    print "(A, I0)", achar(9) // "Function Evaluations: ", rst_%statistics%function_evaluations
    print "(A, I0)", achar(9) // "Matrix Factorizations: ", rst_%statistics%matrix_factorizations
end subroutine

subroutine plot_results()
    type(multiplot) :: plt
    type(plot_2d) :: pltY, pltDYDX, epltY, epltDYDX
    class(legend), pointer :: lgnd, lgnd2
    real(real64), allocatable, dimension(:) :: eAD, eFD, eADdydx, eFDdydx

    call plt%initialize(2, 2, width = 1200, height = 700)
    call pltY%initialize()
    call pltDYDX%initialize()
    call epltY%initialize()
    call epltDYDX%initialize()
    lgnd => pltDYDX%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    lgnd2 => epltDYDX%get_legend()
    call lgnd2%set_is_visible(.true.)
    call lgnd2%set_draw_border(.false.)
    call lgnd2%set_draw_inside_axes(.false.)
    call lgnd2%set_horizontal_position(LEGEND_CENTER)
    call lgnd2%set_vertical_position(LEGEND_BOTTOM)
    call lgnd2%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call plty%set_show_gridlines(.true.)
    call pltDYDX%set_show_gridlines(.true.)
    call epltY%set_show_gridlines(.true.)
    call epltDYDX%set_show_gridlines(.true.)

    call pltY%set_x_axis_title("x")
    call pltY%set_y_axis_title("y")
    call pltY%push(rstAD%solution_points(), rstAD%state(1), name = "AD")
    call pltY%push(rstFD%solution_points(), rstFD%state(1), name = "FD", ls = LINE_DASHED)
    call pltY%push(rst%solution_points(), rst%state(1), name = "Analytical", ls = LINE_DASH_DOT_DOT)
    call plt%set(1, 1, pltY)

    call pltDYDX%set_x_axis_title("x\n")
    call pltDYDX%set_y_axis_title("dy/dx")
    call pltDYDX%push(rstAD%solution_points(), rstAD%state(2), name = "AD")
    call pltDYDX%push(rstFD%solution_points(), rstFD%state(2), name = "FD", ls = LINE_DASHED)
    call pltDYDX%push(rst%solution_points(), rst%state(2), name = "Analytical", ls = LINE_DASH_DOT_DOT)
    call plt%set(2, 1, pltDYDX)

    ! -----
    eAD = rstAD%state(1) - rst%state(1)
    eFD = rstFD%state(1) - rst%state(1)
    eADdydx = rstAD%state(2) - rst%state(2)
    eFDdydx = rstFD%state(2) - rst%state(2)

    call epltY%set_x_axis_title("x")
    call epltY%set_y_axis_title("y_{i} - y_{a}")
    call epltY%push(rst%solution_points(), eAD, name = "AD")
    call epltY%push(rst%solution_points(), eFD, name = "FD", ls = LINE_DASHED)
    call plt%set(1, 2, epltY)

    call epltDYDX%set_x_axis_title("x\n")
    call epltDYDX%set_y_axis_title("dy_{i}/dx - dy_{a}/dx")
    call epltDYDX%push(rst%solution_points(), eADdydx, name = "AD")
    call epltDYDX%push(rst%solution_points(), eFDdydx, name = "FD", ls = LINE_DASHED)
    call plt%set(2, 2, epltDYDX)

    call plt%draw()
end subroutine

end program