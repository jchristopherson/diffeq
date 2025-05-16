program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    implicit none

    ! Initial Conditions & Time Constraints
    real(real64), parameter :: t(2) = [0.0d0, 5.0d1]
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]

    ! Local Variables
    type(runge_kutta_23) :: integrator_1
    type(runge_kutta_45) :: integrator_2
    type(runge_kutta_853) :: integrator_3
    type(rosenbrock) :: integrator_4
    type(bdf) :: integrator_5
    type(adams) :: integrator_6
    type(ode_container) :: mdl
    real(real64), allocatable, dimension(:,:) :: s1, s2, s3, s4, s4a, s5, s6

    ! Define the model
    mdl%fcn => vanderpol

    ! Integrate the model with each integrator
    call integrator_1%solve(mdl, t, ic)
    call integrator_2%solve(mdl, t, ic)
    call integrator_3%solve(mdl, t, ic)
    call integrator_4%solve(mdl, t, ic)
    call integrator_5%solve(mdl, t, ic)
    call integrator_6%solve(mdl, t, ic)

    ! Retrieve the solution from each integrator
    s1 = integrator_1%get_solution()
    s2 = integrator_2%get_solution()
    s3 = integrator_3%get_solution()
    s4 = integrator_4%get_solution()
    s5 = integrator_5%get_solution()
    s6 = integrator_6%get_solution()

    ! Print out the size of each solution
    print "(AI0A)", "RUNGE_KUTTA_23: ", size(s1, 1), " Solution Points"
    print "(AI0A)", "RUNGE_KUTTA_45: ", size(s2, 1), " Solution Points"
    print "(AI0A)", "RUNGE_KUTTA_853: ", size(s3, 1), " Solution Points"
    print "(AI0A)", "ROSENBROCK: ", size(s4, 1), " Solution Points"

    ! Now, implement a PI controller and check its effect.  This will likely
    ! increase the number of steps (loss of efficiency), but if there were
    ! any stability issues, stability will likely improve.  Stability is likely
    ! not relevant on this problem, but it's here for illustration purposes.
    call integrator_4%clear_buffer()
    call integrator_4%set_step_size_control_parameter(0.1d0)
    call integrator_4%solve(mdl, t, ic)
    s4a = integrator_4%get_solution()
    print "(AI0A)", "ROSENBROCK w/ PI Controller: ", size(s4a, 1), " Solution Points"

    ! VODE Integrators
    print "(AI0A)", "BDF: ", size(s5, 1), " Solution Points"
    print "(AI0A)", "ADAMS: ", size(s6, 1), " Solution Points"
end program