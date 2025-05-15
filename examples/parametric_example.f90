program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Initial Conditions & Time Constraints
    real(real64), parameter :: t(2) = [0.0d0, 1.0d0]
    real(real64), parameter :: ic(2) = [0.0d0, 0.0d0]

    ! Local Variables
    type(runge_kutta_853) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define the model
    mdl%fcn => parametric_model

    ! Perform the integration
    call integrator%solve(mdl, t, ic)
    sol = integrator%get_solution()

    ! Plot the solution
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call pd%define_data(sol(:,1), sol(:,2))
    call plt%push(pd)
    call plt%draw()
end program