program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: x(2) = [0.0d0, 5.0d1]
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]

    ! Local Variables
    type(rosenbrock) :: integrator_1
    type(bdf) :: integrator_2
    type(ode_container) :: mdl
    real(real64) :: mu
    real(real64), allocatable, dimension(:,:) :: sol_1, sol_2

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the model
    mdl%fcn => vanderpol_args
    mu = 5.0d0

    ! Compute the solution
    call integrator_1%solve(mdl, x, ic, args = mu)
    sol_1 = integrator_1%get_solution()

    call integrator_2%solve(mdl, x, ic, args = mu)
    sol_2 = integrator_2%get_solution()

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    call lgnd%set_draw_border(.false.)
    
    call pd1%define_data(sol_1(:,1), sol_1(:,2))
    call pd1%set_name("Rosenbrock")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(sol_2(:,1), sol_2(:,2))
    call pd2%set_name("BDF")
    call pd2%set_line_width(3.5)
    call pd2%set_line_style(LINE_DASHED)
    call plt%push(pd2)

    call plt%draw()
end program
