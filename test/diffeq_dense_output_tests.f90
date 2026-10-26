program diffeq_dense_output_tests
    use diffeq
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    implicit none

    class(ode_integrator), allocatable :: solver
    type(ode_container) :: model
    type(ode_results) :: before, after
    character(len=32) :: mode
    integer :: method, direction, output_mode, calls, saved_calls, i, problem
    real(real64) :: x_start, x_end, x_query, y(2), repeated(2), wrong_size(1)
    real(real64), allocatable :: solution(:,:)

    model%fcn => rhs
    call get_command_argument(1, mode)
    if (len_trim(mode) > 0) then
        call error_case(trim(mode))
        stop
    end if

    do method = 1, 9
        select case (method)
        case (1)
            allocate(runge_kutta_23 :: solver)
        case (2)
            allocate(runge_kutta_45 :: solver)
        case (3)
            allocate(runge_kutta_853 :: solver)
        case (4)
            allocate(tsitouras_54 :: solver)
        case (5)
            allocate(rosenbrock :: solver)
        case (6)
            allocate(kennedy_carpenter_4 :: solver)
        case (7)
            allocate(kennedy_carpenter_5 :: solver)
        case (8)
            allocate(bdf :: solver)
        case (9)
            allocate(adams :: solver)
        end select

        call check(.not.solver%get_dense_output_enabled(), "default setting")
        call check(.not.solver%has_dense_output(), "initial availability")
        call solver%set_absolute_tolerance(1.0d-9)
        call solver%set_relative_tolerance(1.0d-9)
        call solver%set_maximum_step_size(0.05d0)
        do direction = -1, 1, 2
            do output_mode = 1, 2
                call solver%clear_buffer()
                call solver%set_dense_output_enabled(output_mode == 1)
                calls = 0
                if (output_mode == 1) then
                    call solver%solve(model, [0.0d0, direction * 0.4d0], &
                        [1.0d0, 1.0d0], args=calls)
                else
                    call solver%solve(model, [0.0d0, direction * 0.2d0, direction * 0.4d0], &
                        [1.0d0, 1.0d0], args=calls)
                end if
                call check(solver%has_dense_output(), "accepted-step availability")
                call solver%get_dense_interval(x_start, x_end)
                call check(direction * (x_end - x_start) > 0.0d0, "interval direction")
                call check(abs(x_end - x_start) <= 0.0500000001d0, "interval is one step")
                before = solver%get_ode_results()
                saved_calls = calls
                do i = 0, 4
                    x_query = x_start + (x_end - x_start) * (i / 4.0d0)
                    if (i == 4) x_query = x_end
                    call solver%evaluate_dense(x_query, y)
                    call check(maxval(abs(y - exact(x_query))) < 3.0d-6, "dense accuracy")
                end do
                call solver%evaluate_dense(x_start, y)
                call solver%evaluate_dense(x_end, repeated)
                solution = solver%get_solution()
                call solver%evaluate_dense(solution(size(solution,1),1), y)
                call check(maxval(abs(y - solution(size(solution,1),2:3))) < 1.0d-12, &
                    "reported endpoint agrees")
                call solver%evaluate_dense(x_end, y)
                call check(all(y == repeated), "repeatability")
                after = solver%get_ode_results()
                call check(calls == saved_calls, "evaluation does not call RHS")
                call check(before%statistics%function_evaluations == &
                    after%statistics%function_evaluations, "evaluation preserves statistics")
                call check(before%count() == after%count(), "evaluation preserves buffer")

                select type (solver)
                class is (multistep_integrator)
                    call solver%interpolate(x_end, repeated)
                    call check(maxval(abs(y - repeated)) < 1.0d-12, "legacy interpolation")
                    call solver%reset_history()
                    call check(.not.solver%has_dense_output(), "history reset invalidates")
                end select
                call solver%clear_buffer()
                call check(.not.solver%has_dense_output(), "buffer reset invalidates")
            end do
        end do

        ! A new solve without preparation must not expose the preceding interpolant.
        call solver%set_dense_output_enabled(.true.)
        call solver%solve(model, [0.0d0, 0.4d0], [1.0d0, 1.0d0], args=calls)
        call check(solver%has_dense_output(), "explicit preparation")
        call solver%set_dense_output_enabled(.false.)
        call check(.not.solver%has_dense_output(), "setting invalidates")
        call solver%clear_buffer()
        call solver%set_maximum_step_size(0.05d0)
        call solver%solve(model, [0.0d0, 0.4d0], [1.0d0, 1.0d0], args=calls)
        if (solver%has_dense_output()) then
            call solver%get_dense_interval(x_start, x_end)
            call check(x_end > 0.4d0, "disabled mode prepares only for overshoot")
        end if
        if (method == 3) then
            before = solver%get_ode_results()
            saved_calls = merge(1, 0, solver%has_dense_output())
            call solver%clear_buffer()
            call solver%set_dense_output_enabled(.true.)
            call solver%solve(model, [0.0d0, 0.4d0], [1.0d0, 1.0d0], args=calls)
            after = solver%get_ode_results()
            call check(before%statistics%successful_steps == after%statistics%successful_steps, &
                "dense preparation preserves step selection")
            call check(after%statistics%function_evaluations - before%statistics%function_evaluations == &
                3 * (after%statistics%successful_steps - saved_calls), "RK853 opt-in cost")
        end if

        if (method >= 5) then
            do problem = 1, 2
                if (problem == 2 .and. method == 9) cycle
                call solver%clear_buffer()
                call solver%set_dense_output_enabled(.true.)
                if (problem == 1) then
                    model%fcn => mass_rhs
                    model%mass_matrix => mass_matrix
                else
                    model%fcn => dae_rhs
                    model%mass_matrix => singular_mass
                end if
                call solver%solve(model, [0.0d0, 0.4d0], [1.0d0, 1.0d0], args=calls)
                call solver%get_dense_interval(x_start, x_end)
                do i = 0, 4
                    x_query = x_start + (x_end - x_start) * (i / 4.0d0)
                    if (i == 4) x_query = x_end
                    call solver%evaluate_dense(x_query, y)
                    call check(maxval(abs(y - exact(x_query))) < 3.0d-6, "mass-matrix dense accuracy")
                end do
            end do
            nullify(model%mass_matrix)
            model%fcn => rhs
        end if

        call solver%clear_buffer()
        call solver%set_dense_output_enabled(.true.)
        call solver%solve(model, [0.0d0, 0.4d0], [1.0d0], args=calls)
        call solver%get_dense_interval(x_start, x_end)
        call solver%evaluate_dense(0.5d0 * (x_start + x_end), wrong_size)
        call check(abs(wrong_size(1) - exp(0.5d0 * (x_start + x_end))) < 3.0d-6, &
            "state-count change")
        deallocate(solver)
    end do
    call check_rk45_coefficients()

contains
    subroutine check_rk45_coefficients()
        type(runge_kutta_45) :: rk45
        type(ode_container) :: polynomial_model
        real(real64) :: y0(2), f0(2), yn(2), fn(2), yerr(2), interpolated(2)
        real(real64), allocatable :: stages(:,:)
        polynomial_model%fcn => polynomial_rhs
        y0 = 0.0d0
        f0 = 0.0d0
        allocate(stages(2, rk45%get_stage_count()))
        call rk45%attempt_step(polynomial_model, 1.0d0, 0.0d0, y0, f0, yn, fn, yerr, stages)
        call rk45%post_step_action(polynomial_model, .true., 0.0d0, 1.0d0, y0, yn, f0, fn, stages)
        call rk45%interpolate(0.5d0, 0.0d0, y0, f0, 1.0d0, yn, fn, interpolated)
        call check(maxval(abs(interpolated - [0.0625d0, 0.125d0])) < 1.0d-12, &
            "RK45 component-specific coefficients")
    end subroutine

    subroutine polynomial_rhs(x, state, derivative, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: derivative(:)
        class(*), intent(inout), optional :: args
        derivative = [4.0d0 * x**3, 8.0d0 * x**3]
    end subroutine

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not.condition) then
            print *, "Dense output test failed: ", label, " method ", method
            error stop 1
        end if
    end subroutine

    function exact(x) result(value)
        real(real64), intent(in) :: x
        real(real64) :: value(2)
        value = [exp(x), exp(-2.0d0 * x)]
    end function

    subroutine rhs(x, state, derivative, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: derivative(:)
        class(*), intent(inout), optional :: args
        if (.not.present(args)) error stop "Missing dense-output user arguments"
        select type (args)
        type is (integer)
            args = args + 1
        class default
            error stop "Unexpected user arguments"
        end select
        derivative(1) = state(1)
        if (size(state) > 1) derivative(2) = -2.0d0 * state(2)
    end subroutine

    subroutine mass_rhs(x, state, derivative, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: derivative(:)
        class(*), intent(inout), optional :: args
        call rhs(x, state, derivative, args)
        derivative = derivative * [2.0d0, 3.0d0]
    end subroutine

    subroutine mass_matrix(x, state, mass, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: mass(:,:)
        class(*), intent(inout), optional :: args
        mass = 0.0d0
        mass(1,1) = 2.0d0
        mass(2,2) = 3.0d0
    end subroutine

    subroutine dae_rhs(x, state, derivative, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: derivative(:)
        class(*), intent(inout), optional :: args
        call rhs(x, state, derivative, args)
        derivative(2) = state(2) - exp(-2.0d0 * x)
    end subroutine

    subroutine singular_mass(x, state, mass, args)
        real(real64), intent(in) :: x, state(:)
        real(real64), intent(out) :: mass(:,:)
        class(*), intent(inout), optional :: args
        mass = 0.0d0
        mass(1,1) = 1.0d0
    end subroutine

    subroutine error_case(name)
        character(len=*), intent(in) :: name
        if (index(name, "multistep") > 0) then
            allocate(bdf :: solver)
        else
            allocate(runge_kutta_853 :: solver)
        end if
        if (index(name, "unavailable") > 0) then
            call solver%evaluate_dense(0.0d0, y)
        else if (index(name, "interval") > 0) then
            call solver%get_dense_interval(x_start, x_end)
        else
            calls = 0
            call solver%set_dense_output_enabled(.true.)
            call solver%solve(model, [0.0d0, -0.4d0], [1.0d0, 1.0d0], args=calls)
            call solver%get_dense_interval(x_start, x_end)
            if (index(name, "size") > 0) then
                call solver%evaluate_dense(x_end, wrong_size)
            else if (index(name, "nan") > 0) then
                call solver%evaluate_dense(ieee_value(0.0d0, ieee_quiet_nan), y)
            else if (index(name, "infinite") > 0) then
                call solver%evaluate_dense(ieee_value(0.0d0, ieee_positive_inf), y)
            else if (index(name, "before") > 0) then
                call solver%evaluate_dense(max(x_start, x_end) + 1.0d0, y)
            else if (index(name, "after") > 0) then
                call solver%evaluate_dense(min(x_start, x_end) - 1.0d0, y)
            else
                error stop "Unknown dense-output error case"
            end if
        end if
    end subroutine
end program
