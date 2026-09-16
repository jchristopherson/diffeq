program diffeq_error_tests
    use diffeq
    implicit none

    character(len=32) :: test_name
    type(runge_kutta_45) :: integrator
    type(ode_container) :: mdl

    call get_command_argument(1, test_name)
    select case (trim(test_name))
    case ("missing_ode")
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("short_grid")
        mdl%fcn => constant_ode
        call integrator%solve(mdl, [0.0d0], [0.0d0])
    case ("step_limit")
        mdl%fcn => constant_ode
        call integrator%set_step_limit(0)
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case ("tolerance_size")
        mdl%fcn => constant_ode
        call integrator%set_absolute_tolerance([1.0d-6, 1.0d-6])
        call integrator%solve(mdl, [0.0d0, 1.0d0], [0.0d0])
    case default
        error stop 2
    end select

contains
    subroutine constant_ode(x, y, dydx, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: dydx(:)
        class(*), intent(inout), optional :: args
        dydx = 1.0d0
    end subroutine
end program