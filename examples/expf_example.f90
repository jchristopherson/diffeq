program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 2500
    real(real64), parameter :: h = 1.0d-2

    ! Local Variables
    type(exponential_fixed_integrator) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts), sol(npts, 3)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define the values of x at which the solution is to be computed
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution
    sol = integrator%solve(mdl, x, [2.0d0, 0.0d0])

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")

    call pd%define_data(sol(:,1), sol(:,2))
    call pd%set_line_width(2.0)
    call plt%push(pd)

    call plt%draw()
end program