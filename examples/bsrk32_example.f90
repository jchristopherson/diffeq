program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: xmax = 2.5d1

    ! Local Variables
    type(bsrk32_integrator) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution
    sol = integrator%solve(mdl, [0.0d0, xmax], [2.0d0, 0.0d0])

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, xmax)

    call pd%define_data(sol(:,1), sol(:,2))
    call pd%set_line_width(2.0)
    call plt%push(pd)

    call plt%draw()
end program