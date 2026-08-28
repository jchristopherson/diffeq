program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Model Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: mass = 2.0d0

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd

    ! Pass in the model parameters
    args%length = length
    args%mass = mass

    ! Define the model
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)

    ! Compute the solution
    call integrator%solve( &
        mdl, &                                  ! model to solve
        [0.0d0, 1.0d1], &                       ! time bounds
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! initial conditions - must satisfy algebraic constraints
        args = args &                           ! arguments to pass to the model
    )
    sol = integrator%get_solution()

! ------------------------------------------------------------------------------
    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()
    call plt%set_use_y2_axis(.true.)
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call y2Axis%set_title("y(t)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_name("x(t)")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(sol(:,1), sol(:,4))
    call pd2%set_name("y(t)")
    call pd2%set_line_width(2.0)
    call pd2%set_draw_against_y2(.true.)
    call plt%push(pd2)

    call plt%draw()
end program