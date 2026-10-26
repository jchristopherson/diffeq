program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Local Variables
    real(real64), parameter :: xmax = 5.0d1
    integer, parameter :: npts = 50
    real(real64) :: x(2), x_dense(npts), ic(2)
    real(real64), allocatable, dimension(:,:) :: sol, dense_sol
    type(kennedy_carpenter_5) :: integrator
    type(ode_container) :: model

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(legend), pointer :: lgnd

    ! Define the x vectors & initial conditions
    ic = [2.0d0, 0.0d0]
    x = [0.0d0, xmax]
    x_dense = linspace(0.0d0, xmax, npts)

    ! Set up the integrator and solve
    model%fcn => vanderpol
    call integrator%solve(model, x, ic)
    sol = integrator%get_solution()

    ! Now, utilize dense output
    call integrator%clear_buffer()  ! need to clear the solution buffer
    call integrator%solve(model, x_dense, ic)   ! simply supply the dense time array
    dense_sol = integrator%get_solution()

! ------------------------------------------------------------------------------
    ! Plot the solutions
    call plt%initialize()
    call plt%set_x_axis_title("x\n")
    call plt%set_y_axis_title("y(x)")
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_line_width(1.0)
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_style(MARKER_FILLED_CIRCLE)
    call pd1%set_marker_scaling(0.6)
    call pd1%set_name("Solution")
    call plt%push(pd1)

    call pd2%define_data(dense_sol(:,1), dense_sol(:,2))
    call pd2%set_draw_line(.false.)
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call pd2%set_name("Dense Output")
    call plt%push(pd2)

    call plt%draw()
end program