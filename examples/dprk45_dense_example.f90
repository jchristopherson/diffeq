program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: xmax = 6.0d1

    ! Local Variables
    type(dprk45_integrator) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: h
    real(real64), allocatable, dimension(:) :: x
    real(real64), allocatable, dimension(:,:) :: std, dense

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the values of x at which the solution is to be computed
    h = xmax / (npts - 1.0d0)
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Define the model
    mdl%fcn => duffing

    ! Compute the solution using standard output and dense output
    std = integrator%solve(mdl, [0.0d0, xmax], [0.0d0, 0.0d0])
    dense = integrator%solve(mdl, x, [0.0d0, 0.0d0])

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call lgnd%set_is_visible(.true.)
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, xmax)

    call pd1%define_data(std(:,1), std(:,2))
    call pd1%set_name("Standard")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(dense(:,1), dense(:,2))
    call pd2%set_name("Dense")
    call pd2%set_line_width(2.0)
    call pd2%set_draw_line(.false.)
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call plt%push(pd2)

    call plt%draw()
end program