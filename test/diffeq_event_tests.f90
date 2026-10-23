program diffeq_event_tests
    use diffeq
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    implicit none

    type(ode_event_solver) :: event_solver
    type(runge_kutta_45) :: integrator
    type(runge_kutta_23) :: integrator23
    type(runge_kutta_853) :: integrator853
    type(tsitouras_54) :: integrator_tsit
    type(ode_container) :: model
    character(len=32) :: test_name
    integer :: callback_count
    real(real64) :: event_x, event_y(1), event_values(2)
    real(real64) :: y1(1), y2(1), x_start, x_end
    real(real64) :: root_a, root_b
    real(real64) :: integration_event_root
    logical :: found
    integer :: event_mode

    call get_command_argument(1, test_name)
    if (len_trim(test_name) > 0) then
        call error_case(trim(test_name))
        error stop "Expected event solver error was not raised"
    end if

    call check(event_solver%get_max_allowed_iterations() == 1000, &
        "default iteration limit")
    call check(event_solver%get_solver_tolerance() == 1.0d-12, &
        "default solver tolerance")
    call check(.not.model%get_is_ode_event_function_defined(), &
        "undefined event callback")
    call check(model%get_event_equation_count() == 0, &
        "default event equation count")

    call event_solver%set_max_allowed_iterations(47)
    call event_solver%set_solver_tolerance(1.0d-9)
    call check(event_solver%get_max_allowed_iterations() == 47, &
        "iteration limit setter")
    call check(event_solver%get_solver_tolerance() == 1.0d-9, &
        "tolerance setter")
    model%event_function => event_values_callback
    call model%set_event_equation_count(2)
    call check(model%get_is_ode_event_function_defined(), &
        "defined event callback")
    call check(model%get_event_equation_count() == 2, &
        "event equation count setter")
    callback_count = 0
    event_mode = 1
    root_a = 0.7d0
    root_b = 0.3d0
    call model%event_function(0.25d0, [0.25d0], event_values, callback_count)
    call check(maxval(abs(event_values - [-0.45d0, -0.05d0])) < 1.0d-14, &
        "event callback evaluation")
    call check(callback_count == 1, "event callback arguments")

    call check(event_solver%is_event_suspected([-1.0d0, 1.0d0], &
        [1.0d0, 2.0d0], 1), "sign-changing bracket")
    call check(event_solver%is_event_suspected([-1.0d0, 0.0d0], &
        [1.0d0, 2.0d0], 2), "zero at bracket endpoint")
    call check(.not.event_solver%is_event_suspected([-1.0d0, 1.0d0], &
        [1.0d0, 2.0d0], 2), "non-crossing bracket")

    model%fcn => unit_rate
    call integrator%set_dense_output_enabled(.true.)
    call integrator%set_absolute_tolerance(1.0d-12)
    call integrator%set_relative_tolerance(1.0d-12)
    call integrator%solve(model, [0.0d0, 1.0d0], [0.0d0])
    call integrator%get_dense_interval(x_start, x_end)
    call integrator%evaluate_dense(x_start, y1)
    call integrator%evaluate_dense(x_end, y2)

    event_mode = 1
    root_a = x_start + 0.7d0 * (x_end - x_start)
    root_b = x_start + 0.3d0 * (x_end - x_start)
    callback_count = 0
    found = event_solver%find_event(model%event_function, integrator, 2, &
        x_start, y1, x_end, y2, event_x, event_y, callback_count)
    call check(found, "forward event found")
    call check(abs(event_x - root_b) < 1.0d-9, "earliest forward event")
    call check(abs(event_y(1) - event_x) < 1.0d-12, "forward event state")
    call check(callback_count > 0, "find event forwards callback arguments")

    event_mode = 3
    root_b = x_start + 0.4d0 * (x_end - x_start)
    found = event_solver%find_event(model%event_function, integrator, 2, &
        x_start, y1, x_end, y2, event_x, event_y)
    call check(found .and. abs(event_x - root_b) < 1.0d-9, &
        "found event selected when another channel has no root")

    call integrator%solve(model, [1.0d0, 0.0d0], [1.0d0])
    call integrator%get_dense_interval(x_start, x_end)
    call integrator%evaluate_dense(x_start, y1)
    call integrator%evaluate_dense(x_end, y2)
    event_mode = 1
    root_a = x_start + 0.7d0 * (x_end - x_start)
    root_b = x_start + 0.3d0 * (x_end - x_start)
    found = event_solver%find_event(model%event_function, integrator, 2, &
        x_start, y1, x_end, y2, event_x, event_y)
    call check(found, "reverse event found")
    call check(abs(event_x - root_b) < 1.0d-9, "earliest reverse event")
    call check(abs(event_y(1) - event_x) < 1.0d-12, "reverse event state")

    event_mode = 2
    found = event_solver%find_event(model%event_function, integrator, 2, &
        x_start, y1, x_end, y2, event_x, event_y)
    call check(.not.found, "no event found")
    call check(ieee_is_nan(event_x) .and. all(ieee_is_nan(event_y)), &
        "no-event outputs are NaNs")

    call check_integrator_event_detection(integrator23, "RK23")
    call check_integrator_event_detection(integrator, "RK45")
    call check_integrator_event_detection(integrator853, "RK853")
    call check_integrator_event_detection(integrator_tsit, "Tsitouras")

contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not.condition) then
            print *, "Event solver test failed: ", label
            error stop 1
        end if
    end subroutine

    subroutine unit_rate(x, y, dydx, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: dydx(:)
        class(*), intent(inout), optional :: args
        dydx = 1.0d0
    end subroutine

    subroutine check_integrator_event_detection(solver, solver_name)
        class(single_step_integrator), intent(inout) :: solver
        character(len=*), intent(in) :: solver_name
        real(real64), allocatable :: solution(:,:)
        real(real64) :: root

        model%fcn => unit_rate
        model%event_function => single_event_callback
        call model%set_event_equation_count(1)

        ! The default two-point output grid does not otherwise request dense
        ! output; event localization must still use the accepted step interpolant.
        call solver%clear_buffer()
        call solver%set_dense_output_enabled(.false.)
        integration_event_root = 1.5d-6
        call solver%solve(model, [0.0d0, 1.0d0], [0.0d0])
        solution = solver%get_solution()
        call check(size(solution, 1) >= 2, &
            solver_name // " two-point event output size")
        call check(abs(solution(size(solution,1),1) - &
            integration_event_root) < 1.0d-10, &
            solver_name // " two-point event location")
        call check(abs(solution(size(solution,1),2) - &
            integration_event_root) < 1.0d-10, &
            solver_name // " two-point event state")
        call check(all(solution(:,1) <= integration_event_root + 1.0d-10), &
            solver_name // " two-point output stops at event")

        ! A requested output point before the event may share its accepted
        ! step, but must remain in the solution returned at termination.
        call solver%clear_buffer()
        call solver%set_dense_output_enabled(.false.)
        integration_event_root = 1.5d-6
        call solver%solve(model, [0.0d0, 1.0d-6, 1.0d0], [0.0d0])
        solution = solver%get_solution()
        call check(size(solution, 1) == 3, &
            solver_name // " keeps pre-event output")
        call check(abs(solution(2,1) - 1.0d-6) < 1.0d-12, &
            solver_name // " requested output before event")
        call check(abs(solution(2,2) - 1.0d-6) < 1.0d-10, &
            solver_name // " requested output state before event")
        call check(abs(solution(3,1) - integration_event_root) < 1.0d-10, &
            solver_name // " event follows pre-event output")
        call check(abs(solution(3,2) - integration_event_root) < 1.0d-10, &
            solver_name // " event state after pre-event output")

        ! Repeat in reverse integration order to verify ordering and
        ! localization are independent of the sign of the step.
        call solver%clear_buffer()
        call solver%set_dense_output_enabled(.false.)
        root = 1.0d0 - 1.5d-6
        integration_event_root = root
        call solver%solve(model, [1.0d0, 1.0d0 - 1.0d-6, 0.0d0], [1.0d0])
        solution = solver%get_solution()
        call check(size(solution, 1) == 3, &
            solver_name // " reverse keeps pre-event output")
        call check(abs(solution(2,1) - (1.0d0 - 1.0d-6)) < 1.0d-12, &
            solver_name // " reverse requested output before event")
        call check(abs(solution(3,1) - root) < 1.0d-10, &
            solver_name // " reverse event location")
        call check(abs(solution(3,2) - root) < 1.0d-10, &
            solver_name // " reverse event state")
    end subroutine

    subroutine single_event_callback(x, y, values, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: values(:)
        class(*), intent(inout), optional :: args
        values = x - integration_event_root
    end subroutine

    subroutine event_values_callback(x, y, values, args)
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: values(:)
        class(*), intent(inout), optional :: args
        if (present(args)) then
            select type (args)
            type is (integer)
                args = args + 1
            class default
                error stop "Unexpected event callback argument"
            end select
        end if
        select case (event_mode)
        case (1)
            values = [x - root_a, x - root_b]
        case (2)
            values = [1.0d0, 2.0d0]
        case (3)
            values = [1.0d0, x - root_b]
        case default
            error stop "Unexpected event mode"
        end select
    end subroutine

    subroutine error_case(name)
        character(len=*), intent(in) :: name
        real(real64) :: x
        real(real64) :: y(1)
        real(real64) :: wrong_y(2)
        logical :: status

        wrong_y = 0.0d0
        select case (name)
        case ("bracket_size")
            status = event_solver%is_event_suspected([1.0d0], [1.0d0, 2.0d0], 1)
        case ("bracket_index")
            status = event_solver%is_event_suspected([1.0d0], [1.0d0], 2)
        case default
            model%event_function => event_values_callback
            select case (name)
            case ("neqn")
                status = event_solver%find_event(model%event_function, &
                    integrator, 0, 0.0d0, [0.0d0], 1.0d0, [1.0d0], x, y)
            case ("equal_bounds")
                status = event_solver%find_event(model%event_function, &
                    integrator, 1, 0.0d0, [0.0d0], 0.0d0, [0.0d0], x, y)
            case ("y2_size")
                status = event_solver%find_event(model%event_function, &
                    integrator, 1, 0.0d0, [0.0d0], 1.0d0, [1.0d0, 2.0d0], &
                    x, y)
            case ("y_size")
                status = event_solver%find_event(model%event_function, &
                    integrator, 1, 0.0d0, [0.0d0], 1.0d0, [1.0d0], x, &
                    wrong_y)
            case default
                error stop "Unknown event solver error case"
            end select
        end select
    end subroutine
end program
