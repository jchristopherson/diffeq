program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 25000
    real(real64), parameter :: h = 1.0d-3

    ! Local Variables
    type(adams_fixed_integerator) :: integrator
    type(rk4_fixed_integrator) :: comparison
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts), sol(npts, 3), csol(npts, 3)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the values of x at which the solution is to be computed
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution with both integrators
    sol = integrator%solve(mdl, x, [2.0d0, 0.0d0])
    csol = comparison%solve(mdl, x, [2.0d0, 0.0d0])

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_name("Adams")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(csol(:,1), csol(:,2))
    call pd2%set_name("Exponential")
    call pd2%set_line_width(4.0)
    call pd2%set_line_style(LINE_DASHED)
    call plt%push(pd2)

    call plt%draw()
end program