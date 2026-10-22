program example
    use iso_fortran_env
    use diffeq
    use fplot_core
    implicit none

    ! Local Variables
    real(real64), parameter :: tmax = 2.5d2
    integer, parameter :: npts = 50
    real(real64) :: t(2), t_dense(npts), ic(2)
    real(real64), allocatable, dimension(:,:) :: sol, dense_sol
    type(kennedy_carpenter_5) :: integrator
    type(ode_container) :: model

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(legend), pointer :: lgnd

    ! Define the time vectors & initial conditions
    ic = [2.0d0, 0.0d0]
    t = [0.0d0, tmax]
    t_dense = linspace(0.0d0, tmax, npts)

    ! Set up the integrator and solve
    model%fcn => eqn
    call integrator%solve(model, t, ic)
    sol = integrator%get_solution()

    ! Now, utilize dense output
    call integrator%clear_buffer()  ! need to clear the solution buffer
    call integrator%solve(model, t_dense, ic)   ! simply supply the dense time array
    dense_sol = integrator%get_solution()

! ------------------------------------------------------------------------------
    ! Plot the solutions
    call plt%initialize()
    call plt%set_x_axis_title("t")
    call plt%set_y_axis_title("x(t)")
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

contains
    subroutine eqn(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args

        real(real64), parameter :: mu = 5.0d1

        dxdt(1) = x(2)
        dxdt(2) = mu * (1.0d0 - x(1)**2) * x(2) - x(1)
    end subroutine
end program