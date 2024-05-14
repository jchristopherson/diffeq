program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution
    call integrator%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
    sol = integrator%get_solution()

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call y2Axis%set_title("y'(x)")
    call plt%set_use_y2_axis(.true.)
    call lgnd%set_is_visible(.true.)
    
    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_name("y(x)")
    call plt%push(pd1)

    call pd2%define_data(sol(:,1), sol(:,3))
    call pd2%set_draw_against_y2(.true.)
    call pd2%set_name("y'(x)")
    call plt%push(pd2)

    call plt%draw()
end program