program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: xmax = 2.5d1

    ! Local Variables
    type(sdirk4_integrator) :: integrator
    type(dprk45_integrator) :: ref_integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ref(:,:)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution
    call integrator%set_use_pi_controller(.false.)
    sol = integrator%solve(mdl, [0.0d0, xmax], [2.0d0, 0.0d0])
    ref = ref_integrator%solve(mdl, [0.0d0, xmax], [2.0d0, 0.0d0])

    ! Compare the number of solution points for each integrator
    print *, "SDIRK: ", size(sol, 1)
    print *, "DPRK45: ", size(ref, 1)

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, xmax)
    call lgnd%set_is_visible(.true.)

    call pd%define_data(sol(:,1), sol(:,2))
    call pd%set_name("SDIRK4")
    call pd%set_line_width(2.0)
    call plt%push(pd)

    call pd%define_data(ref(:,1), ref(:,2))
    call pd%set_name("DPRK45")
    call pd%set_line_width(2.0)
    call pd%set_line_style(LINE_DASHED)
    call plt%push(pd)

    call plt%draw()
end program