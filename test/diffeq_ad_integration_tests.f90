program diffeq_ad_integration_tests
    !! Container dispatch, solver accuracy, and exact RHS accounting for AD.
    use diffeq
    use diffeq_auto_diff
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none

    type callback_counts
        !! User-owned callback data; counters do not affect the equations.
        integer :: rhs = 0
        integer :: jac = 0
        integer :: state_seeds = 0
        integer :: time_seeds = 0
        logical :: dae = .false.
        logical :: variable_mass = .false.
    end type

    character(40) :: mode
    call get_command_argument(1, mode)
    if (len_trim(mode) > 0) then
        call error_case(trim(mode))
        stop
    end if
    call test_container()
    call test_solvers()
    call test_dae()
    print *, "AD integration tests passed."

contains

    subroutine check(condition, label)
        !! Fail with the expectation that was not met.
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not.condition) then
            print *, "FAIL: ", label
            error stop 1
        end if
    end subroutine

    subroutine near(actual, expected, tolerance, label)
        !! Check a finite array result against an absolute accuracy threshold.
        real(real64), intent(in) :: actual(:), expected(:), tolerance
        character(*), intent(in) :: label
        call check(size(actual) == size(expected), label // " size")
        call check(all(ieee_is_finite(actual)), label // " finite")
        if (maxval(abs(actual - expected)) > tolerance) then
            print *, "FAIL: ", label, " actual: ", actual, " expected: ", expected
            error stop 1
        end if
    end subroutine

    subroutine nonlinear_ad(x, y, f, args)
        !! Coupled nonlinear expression with explicit time dependence.
        type(real_ad), intent(in) :: x, y(:)
        type(real_ad), intent(out) :: f(:)
        class(*), intent(inout), optional :: args
        integer :: i
        if (present(args)) then
            select type (args)
            type is (callback_counts)
                args%rhs = args%rhs + 1
                call check(count(y%dx /= 0.0d0) <= 1, "only one state seeded")
                if (x%dx /= 0.0d0) then
                    call check(all(y%dx == 0.0d0), "time differentiation holds states fixed")
                    args%time_seeds = args%time_seeds + 1
                end if
            class default
                error stop "Unexpected nonlinear callback data"
            end select
        end if
        do i = 1, size(y)
            f(i) = x * y(i)**2 + sin(y(1))
        end do
    end subroutine

    subroutine test_container()
        !! Verify dispatch, analytic precedence, seed reset, and workspace resizing.
        type(ode_container) :: model
        type(callback_counts) :: counts
        real(real64), allocatable :: y(:), f(:), dt(:), jac(:,:), exact(:,:)
        integer :: n, i, j, calls
        real(real64) :: small_jac(1,1), pair_jac(2,2), pair_rhs(2), pair_dt(2)
        real(real64), parameter :: x = 0.3d0

        call check(.not.model%get_is_ode_defined(), "empty model")
        call check(model%get_jacobian_function_evaluations() == 0, "initial Jacobian count")
        model%fcn_ad => nonlinear_ad
        call check(model%get_is_ode_defined(), "AD-only model is defined")
        call check(.not.model%use_finite_difference_jacobian(), "AD is not finite differences")
        call model%set_finite_difference_step(0.125d0)
        do n = 1, 3
            allocate(y(n), f(n), dt(n), jac(n,n), exact(n,n))
            y = [(0.2d0 * i, i = 1, n)]
            counts = callback_counts()
            call model%compute_rhs(x, y, f, counts)
            call near(f, x * y**2 + sin(y(1)), 1.0d-14, "AD values")
            call check(counts%rhs == 1, "one RHS evaluation")
            exact = 0.0d0
            exact(:,1) = cos(y(1))
            do i = 1, n
                exact(i,i) = exact(i,i) + 2.0d0 * x * y(i)
            end do
            do j = 1, 2
                calls = counts%rhs
                call model%compute_jacobian(x, y, jac, counts)
                call near(reshape(jac, [n*n]), reshape(exact, [n*n]), 1.0d-14, "AD Jacobian")
                call check(counts%rhs - calls == n, "one callback per column")
                call check(model%get_jacobian_function_evaluations() == n, "reported AD calls")
                calls = counts%rhs
                call model%compute_time_derivative(x, y, f, dt, counts)
                call near(dt, y**2, 1.0d-14, "AD time derivative")
                call check(counts%rhs == calls + 1, "one callback for time derivative")
                call model%compute_rhs(x, y, f, counts)
                call near(f, x * y**2 + sin(y(1)), 1.0d-14, "values after seeded evaluation")
            end do
            call model%compute_rhs(x, y, f)
            call model%compute_jacobian(x, y, jac)
            call model%compute_time_derivative(x, y, f, dt)
            deallocate(y, f, dt, jac, exact)
        end do
        ! Shrink the same workspace and assign a constant output to clear its seed.
        call model%compute_jacobian(x, [0.2d0], small_jac, counts)
        call near([small_jac(1,1)], [0.12d0 + cos(0.2d0)], 1.0d-14, "shrunk AD workspace")
        model%fcn_ad => constant_ad
        call model%compute_jacobian(x, [0.2d0], small_jac)
        call near([small_jac(1,1)], [0.0d0], 0.0d0, "constant output clears derivative")

        do i = 1, 4
            call configure(model, i)
            counts = callback_counts()
            call model%compute_jacobian(x, [1.0d0, 1.0d0], pair_jac, counts)
            call near(reshape(pair_jac, [4]), [-2.0d0, 0.0d0, 0.0d0, -3.0d0], 1.0d-7, "selected Jacobian")
            select case (i)
            case (1)
                call check(model%use_finite_difference_jacobian(), "real RHS selects FD")
                call check(counts%rhs == 3, "FD callback cost")
            case (2)
                call check(counts%rhs == 2, "AD callback cost")
            case default
                call check(counts%rhs == 0 .and. counts%jac == 1, "explicit Jacobian precedence")
            end select
            call check(model%get_jacobian_function_evaluations() == counts%rhs, "reported selected cost")
            call model%compute_rhs(x, [1.0d0, 1.0d0], pair_rhs, counts)
            call model%compute_time_derivative(x, [1.0d0, 1.0d0], pair_rhs, pair_dt, counts)
            call near(pair_dt, [-exp(-x), -2.0d0 * exp(-x)], 1.0d-7, "selected time derivative")
        end do
        nullify(model%jacobian)
        call model%compute_jacobian(x, [1.0d0, 1.0d0], pair_jac, counts)
        call check(model%get_jacobian_function_evaluations() == 2, "return to AD resets count")
        nullify(model%fcn_ad)
        model%fcn => model_real
        call model%compute_jacobian(x, [1.0d0, 1.0d0], pair_jac, counts)
        call check(model%get_jacobian_function_evaluations() == 3, "switch to FD resets count")
        model%jacobian => model_jac
        call model%compute_jacobian(x, [1.0d0, 1.0d0], pair_jac, counts)
        call check(model%get_jacobian_function_evaluations() == 0, "return to explicit resets count")
    end subroutine

    subroutine constant_ad(x, y, f, args)
        !! Assignment of a constant clears both output derivative components.
        type(real_ad), intent(in) :: x, y(:)
        type(real_ad), intent(out) :: f(:)
        class(*), intent(inout), optional :: args
        f = 2.0d0
    end subroutine

    subroutine make_solver(method, solver)
        !! Construct each supported integrator through its existing base API.
        integer, intent(in) :: method
        class(ode_integrator), allocatable, intent(out) :: solver
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
        case default
            error stop "Unknown method"
        end select
        call solver%set_absolute_tolerance(1.0d-9)
        call solver%set_relative_tolerance(1.0d-9)
        call solver%set_maximum_step_size(0.02d0)
        call solver%set_dense_output_enabled(.true.)
    end subroutine

    subroutine configure(model, representation)
        !! Select FD, AD, real plus analytic, or AD plus analytic callbacks.
        type(ode_container), intent(out) :: model
        integer, intent(in) :: representation
        if (representation == 2 .or. representation == 4) then
            model%fcn_ad => model_ad
        else
            model%fcn => model_real
        end if
        if (representation >= 3) model%jacobian => model_jac
    end subroutine

    subroutine model_real(x, y, f, args)
        !! Manufactured ODE/DAE with known exponential solution.
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: f(:)
        class(*), intent(inout), optional :: args
        select type (args)
        type is (callback_counts)
            args%rhs = args%rhs + 1
            if (args%dae) then
                f(1) = -y(1)
                if (args%variable_mass) f(1) = -(1.0d0 + x + y(1)) * y(1)
                f(2) = y(1)**2 + exp(-2.0d0 * x) - 2.0d0 * y(2)
            else
                f(1) = -2.0d0 * y(1) + exp(-x)
                f(2) = -3.0d0 * y(2) + 2.0d0 * exp(-x)
            end if
        class default
            error stop "Expected callback_counts"
        end select
    end subroutine

    subroutine model_ad(x, y, f, args)
        !! AD representation, counting every seed evaluation independently.
        type(real_ad), intent(in) :: x, y(:)
        type(real_ad), intent(out) :: f(:)
        class(*), intent(inout), optional :: args
        select type (args)
        type is (callback_counts)
            args%rhs = args%rhs + 1
            if (any(y%dx /= 0.0d0)) args%state_seeds = args%state_seeds + 1
            if (x%dx /= 0.0d0) args%time_seeds = args%time_seeds + 1
            if (args%dae) then
                f(1) = -y(1)
                if (args%variable_mass) f(1) = -(1.0d0 + x + y(1)) * y(1)
                f(2) = y(1)**2 + exp(-2.0d0 * x) - 2.0d0 * y(2)
            else
                f(1) = -2.0d0 * y(1) + exp(-x)
                f(2) = -3.0d0 * y(2) + 2.0d0 * exp(-x)
            end if
        class default
            error stop "Expected callback_counts"
        end select
    end subroutine

    subroutine model_jac(x, y, jac, args)
        !! Explicit Jacobian must override the AD or FD differentiation path.
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: jac(:,:)
        class(*), intent(inout), optional :: args
        select type (args)
        type is (callback_counts)
            args%jac = args%jac + 1
            jac = 0.0d0
            if (args%dae) then
                jac(1,1) = -1.0d0
                if (args%variable_mass) jac(1,1) = -(1.0d0 + x + 2.0d0 * y(1))
                jac(2,1) = 2.0d0 * y(1)
                jac(2,2) = -2.0d0
            else
                jac(1,1) = -2.0d0
                jac(2,2) = -3.0d0
            end if
        class default
            error stop "Expected callback_counts"
        end select
    end subroutine

    subroutine mass(x, y, matrix, args)
        !! Singular mass matrix, optionally dependent on both time and state.
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: matrix(:,:)
        class(*), intent(inout), optional :: args
        select type (args)
        type is (callback_counts)
            matrix = 0.0d0
            matrix(1,1) = 1.0d0
            if (args%variable_mass) matrix(1,1) = 1.0d0 + x + y(1)
        class default
            error stop "Expected callback_counts"
        end select
    end subroutine

    subroutine check_statistics(result, counts, representation)
        !! Compare solver statistics to observed callbacks, not estimated costs.
        type(ode_results), intent(in) :: result
        type(callback_counts), intent(in) :: counts
        integer, intent(in) :: representation
        call check(result%statistics%function_evaluations == counts%rhs, "exact total RHS statistics")
        if (representation >= 3) then
            call check(result%statistics%jacobian_evaluations == counts%jac, "analytic Jacobian statistics")
            call check(counts%state_seeds == 0, "explicit Jacobian suppresses AD columns")
        else if (representation == 2) then
            call check(counts%state_seeds == 2 * result%statistics%jacobian_evaluations, &
                "two AD callbacks per Jacobian")
        end if
    end subroutine

    subroutine test_solvers()
        !! Cover every integrator, all callback modes, dense evaluation and events.
        class(ode_integrator), allocatable :: solver
        type(ode_container) :: model
        type(callback_counts) :: counts
        type(ode_results) :: result
        real(real64), allocatable :: solution(:,:)
        real(real64) :: left, right, query, dense(2)
        integer :: method, representation, saved

        do method = 1, 9
            do representation = 1, 4
                print *, "ODE method/mode: ", method, representation
                call make_solver(method, solver)
                call configure(model, representation)
                counts = callback_counts()
                call solver%solve(model, [0.0d0, 0.2d0], [1.0d0, 1.0d0], args=counts)
                result = solver%get_ode_results()
                solution = solver%get_solution()
                call near(solution(size(solution,1),2:3), [exp(-0.2d0), exp(-0.2d0)], 2.0d-6, "ODE accuracy")
                call check_statistics(result, counts, representation)
                if (method == 5 .and. (representation == 2 .or. representation == 4)) &
                    call check(counts%time_seeds > 0, "Rosenbrock AD time derivative")
                call solver%get_dense_interval(left, right)
                query = (left + right) / 2.0d0
                saved = counts%rhs
                call solver%evaluate_dense(query, dense)
                call near(dense, [exp(-query), exp(-query)], 3.0d-6, "AD dense accuracy")
                call check(counts%rhs == saved, "dense output makes no RHS calls")
                deallocate(solver)
            end do

            call make_solver(method, solver)
            call configure(model, 2)
            model%event_function => event
            call model%set_event_equation_count(1)
            counts = callback_counts()
            call solver%solve(model, [0.0d0, 0.4d0], [1.0d0, 1.0d0], args=counts)
            result = solver%get_ode_results()
            solution = solver%get_solution()
            call near([solution(size(solution,1),1)], [0.15d0], 5.0d-6, "AD event time")
            call check_statistics(result, counts, 2)
            deallocate(solver)
        end do
    end subroutine

    subroutine event(x, y, values, args)
        !! Terminal crossing at x = 0.15 along the manufactured ODE solution.
        real(real64), intent(in) :: x, y(:)
        real(real64), intent(out) :: values(:)
        class(*), intent(inout), optional :: args
        values(1) = y(1) - exp(-0.15d0)
    end subroutine

    subroutine test_dae()
        !! Exercise constant/variable-mass DAE Jacobians and initialization counts.
        class(ode_integrator), allocatable :: solver
        type(ode_container) :: model
        type(callback_counts) :: counts
        type(ode_results) :: result
        real(real64), allocatable :: solution(:,:)
        integer :: method, representation, dependent

        do method = 5, 8
            do dependent = 0, 1
                do representation = 1, 4
                    print *, "DAE method/mass/mode: ", method, dependent, representation
                    call make_solver(method, solver)
                    call configure(model, representation)
                    counts = callback_counts(dae=.true., variable_mass=dependent == 1)
                    model%mass_matrix => mass
                    call model%set_is_mass_matrix_dependent(dependent == 1)
                    call solver%solve(model, [0.0d0, 0.1d0], [1.0d0, 1.0d0], args=counts)
                    result = solver%get_ode_results()
                    solution = solver%get_solution()
                    call near(solution(size(solution,1),2:3), [exp(-0.1d0), exp(-0.2d0)], 3.0d-6, "DAE accuracy")
                    call check_statistics(result, counts, representation)
                    if ((method == 6 .or. method == 7) .and. &
                        (representation == 2 .or. representation == 4)) &
                        call check(counts%time_seeds > 0, "KC DAE AD time derivative")
                    deallocate(solver)
                end do
            end do
        end do
    end subroutine

    subroutine error_case(name)
        !! Negative tests use repository-standard error-stop codes and CTest WILL_FAIL.
        character(*), intent(in) :: name
        type(ode_container) :: model
        type(ggl_container) :: ggl
        type(runge_kutta_45) :: solver
        type(callback_counts) :: counts
        real(real64) :: f(2), dt(2), jac(2,2)
        logical :: defined

        model%fcn_ad => nonlinear_ad
        select case (name)
        case ("both_rhs", "both_jacobian", "both_solve", "both_inquiry", "both_time")
            model%fcn => model_real
            select case (name)
            case ("both_rhs")
                call model%compute_rhs(0.0d0, [1.0d0, 1.0d0], f, counts)
            case ("both_jacobian")
                model%jacobian => model_jac
                call model%compute_jacobian(0.0d0, [1.0d0, 1.0d0], jac, counts)
            case ("both_solve")
                call solver%solve(model, [0.0d0, 0.1d0], [1.0d0, 1.0d0], args=counts)
            case ("both_inquiry")
                defined = model%get_is_ode_defined()
            case ("both_time")
                call model%compute_time_derivative(0.0d0, [1.0d0, 1.0d0], [0.0d0, 0.0d0], dt)
            end select
        case ("rhs_size")
            call model%compute_rhs(0.0d0, [1.0d0], f)
        case ("jacobian_size")
            call model%compute_jacobian(0.0d0, [1.0d0], jac)
        case ("time_size")
            call model%compute_time_derivative(0.0d0, [1.0d0], [0.0d0], dt)
        case ("missing_rhs")
            nullify(model%fcn_ad)
            call model%compute_rhs(0.0d0, [1.0d0, 1.0d0], f)
        case ("ggl_ad", "ggl_ad_jacobian")
            ggl%fcn_ad => nonlinear_ad
            if (name == "ggl_ad") then
                defined = ggl%get_is_ode_defined()
            else
                ggl%jacobian => model_jac
                call ggl%compute_jacobian(0.0d0, [1.0d0, 1.0d0], jac, counts)
            end if
        case default
            error stop "Unknown AD error case"
        end select
    end subroutine
end program
