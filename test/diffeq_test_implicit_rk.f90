! This file is part of diffeq.
! 
! diffeq is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! diffeq is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with diffeq. If not, see <https://www.gnu.org/licenses/>.

module diffeq_test_implicit_rk
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

    type jacobian_counter
        integer :: calls = 0
        integer :: mass_calls = 0
    end type

    type, extends(ode_container) :: mass_dispatch_model
        logical :: mass_enabled = .false.
        integer :: mass_calls = 0
    contains
        procedure :: has_mass_matrix => dispatch_has_mass_matrix
        procedure :: compute_mass_matrix => dispatch_mass_matrix
    end type

contains
! ------------------------------------------------------------------------------
subroutine configure_ggl_test_model(mdl)
    type(ggl_container), intent(out) :: mdl

    call mdl%initialize(2, 1, alpha = 2.0d0, beta = 3.0d0)
    mdl%fcn => ggl_pendulum_force
    mdl%mass_matrix => ggl_pendulum_mass
    mdl%constraints => ggl_pendulum_constraint
    mdl%constraint_jacobian => ggl_pendulum_gradient
    mdl%constraint_curvature => ggl_pendulum_curvature
end subroutine

subroutine ggl_pendulum_force(x, y, force, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: force
    class(*), intent(inout), optional :: args

    select type (args)
    type is (cartesian_pendulum_properties)
        force = [0.0d0, -args%mass * 9.81d0]
    class default
        error stop "Missing pendulum properties"
    end select
end subroutine

subroutine ggl_pendulum_mass(x, q, mass, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:,:) :: mass
    class(*), intent(inout), optional :: args

    select type (args)
    type is (cartesian_pendulum_properties)
        mass = 0.0d0
        mass(1,1) = args%mass
        mass(2,2) = args%mass
    class default
        error stop "Missing pendulum properties"
    end select
end subroutine

subroutine ggl_pendulum_variable_mass(x, q, mass, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:,:) :: mass
    class(*), intent(inout), optional :: args

    call ggl_pendulum_mass(x, q, mass, args)
    mass = mass * (1.0d0 + x)
end subroutine

subroutine ggl_pendulum_constraint(x, q, g, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:) :: g
    class(*), intent(inout), optional :: args

    select type (args)
    type is (cartesian_pendulum_properties)
        g(1) = dot_product(q, q) - args%length**2
    class default
        error stop "Missing pendulum properties"
    end select
end subroutine

subroutine ggl_pendulum_gradient(x, q, gradient, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:,:) :: gradient
    class(*), intent(inout), optional :: args

    gradient(1,:) = 2.0d0 * q
end subroutine

subroutine ggl_pendulum_curvature(x, q, direction, v, curvature, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: q, direction, v
    real(real64), intent(out), dimension(:) :: curvature
    class(*), intent(inout), optional :: args

    curvature(1) = 2.0d0 * dot_product(direction, v)
end subroutine

function test_ggl_container() result(rst)
    logical :: rst
    type(ggl_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), dimension(6) :: state, rhs
    real(real64), dimension(6,6) :: mass, jac

    rst = mdl%has_mass_matrix() .and. .not.mdl%get_is_ode_defined()
    rst = rst .and. mdl%get_state_count() == 0
    call mdl%initialize(2, 1)
    mdl%fcn => ggl_pendulum_force
    mdl%constraints => ggl_pendulum_constraint
    mdl%constraint_jacobian => ggl_pendulum_gradient
    rst = rst .and. .not.mdl%get_is_ode_defined()
    mdl%constraint_curvature => ggl_pendulum_curvature
    rst = rst .and. mdl%get_is_ode_defined() .and. mdl%get_state_count() == 6
    args%length = 1.5d0
    args%mass = 2.0d0
    state = [args%length, 0.0d0, 0.0d0, 0.5d0, 0.0d0, 0.0d0]
    call mdl%compute_mass_matrix(0.0d0, state, mass, args)
    rst = rst .and. all(abs(mass(5:6,:)) < epsilon(1.0d0))
    rst = rst .and. all(abs(mass(:,5:6)) < epsilon(1.0d0))
    rst = rst .and. mass(1,1) == 1.0d0 .and. mass(3,3) == 1.0d0
    call mdl%initialize_multipliers(0.0d0, state, args)
    rst = rst .and. abs(state(5) - 0.25d0 / (2.0d0 * args%length**2)) < 1.0d-12
    rst = rst .and. mdl%validate_initial_conditions(0.0d0, state, 1.0d-12, args)
    mdl%mass_matrix => ggl_pendulum_mass
    call mdl%initialize_multipliers(0.0d0, state, args)
    call mdl%compute_mass_matrix(0.0d0, state, mass, args)
    rst = rst .and. mass(3,3) == args%mass .and. mass(4,4) == args%mass
    rst = rst .and. abs(state(5) - args%mass * 0.25d0 / &
        (2.0d0 * args%length**2)) < 1.0d-12
    call mdl%compute_rhs(0.0d0, state, rhs, args)
    rst = rst .and. maxval(abs(rhs(5:6))) < 1.0d-12
    rst = rst .and. all(state(1:4) == [args%length, 0.0d0, 0.0d0, 0.5d0])
    call mdl%compute_jacobian(0.0d0, state, jac, args)
    rst = rst .and. abs(jac(5,6) - 4.0d0 * args%length**2) < 1.0d-6
    rst = rst .and. abs(jac(6,5) + 4.0d0 * args%length**2 / args%mass) < 1.0d-6
    state(3) = 0.01d0
    rst = rst .and. .not.mdl%validate_initial_conditions(0.0d0, state, 1.0d-6, args)
    state(3) = 0.0d0
    state(1) = state(1) + 0.01d0
    rst = rst .and. .not.mdl%validate_initial_conditions(0.0d0, state, 1.0d-6, args)
    if (.not.rst) print *, "TEST FAILED: test_ggl_container"
end function

function check_ggl_solver(integrator) result(rst)
    class(ode_integrator), intent(inout) :: integrator
    logical :: rst
    type(ggl_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), dimension(6) :: state, rhs
    real(real64), dimension(5) :: times
    real(real64) :: initial_position, initial_velocity, alpha, beta, position, velocity
    real(real64), allocatable, dimension(:,:) :: sol
    integer :: scenario, point

    rst = .true.
    args%length = 1.5d0
    args%mass = 2.0d0
    times = [0.0d0, 0.025d0, 0.05d0, 0.075d0, 0.1d0]
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%set_maximum_step_size(0.01d0)
    do scenario = 1, 4
        call configure_ggl_test_model(mdl)
        alpha = 2.0d0
        beta = 3.0d0
        state = [args%length, 0.0d0, 0.0d0, 0.2d0, 0.0d0, 0.0d0]
        if (scenario > 1) then
            state(1) = state(1) + 0.01d0
            state(3) = 0.01d0
        end if
        if (scenario == 3) then
            alpha = 0.0d0
            beta = 0.0d0
        end if
        if (scenario == 4) then
            mdl%mass_matrix => ggl_pendulum_variable_mass
            call mdl%set_is_mass_matrix_dependent(.true.)
        end if
        call mdl%set_stabilization_parameters(alpha, beta)
        initial_position = dot_product(state(1:2), state(1:2)) - args%length**2
        initial_velocity = 2.0d0 * dot_product(state(1:2), state(3:4))
        call mdl%initialize_multipliers(times(1), state, args)
        call mdl%compute_rhs(times(1), state, rhs, args)
        rst = rst .and. maxval(abs(rhs(5:6))) < 1.0d-12
        call integrator%clear_buffer()
        call integrator%solve(mdl, times, state, args)
        sol = integrator%get_solution()
        do point = 1, size(sol,1)
            position = dot_product(sol(point,2:3), sol(point,2:3)) - args%length**2
            velocity = 2.0d0 * dot_product(sol(point,2:3), sol(point,4:5))
            rst = rst .and. abs(position - initial_position * exp(-alpha * sol(point,1))) < 2.0d-6
            rst = rst .and. abs(velocity - initial_velocity * exp(-beta * sol(point,1))) < 2.0d-6
        end do
    end do
end function

function test_ggl_integration() result(rst)
    logical :: rst
    logical, dimension(4) :: results
    type(bdf) :: bdf_integrator
    type(rosenbrock) :: rosenbrock_integrator
    type(kennedy_carpenter_4) :: kc4_integrator
    type(kennedy_carpenter_5) :: kc5_integrator

    results(1) = check_ggl_solver(bdf_integrator)
    results(2) = check_ggl_solver(rosenbrock_integrator)
    results(3) = check_ggl_solver(kc4_integrator)
    results(4) = check_ggl_solver(kc5_integrator)
    rst = all(results)
    if (.not.rst) print *, "TEST FAILED: test_ggl_integration", results
end function

! ------------------------------------------------------------------------------
function dispatch_has_mass_matrix(this) result(rst)
    class(mass_dispatch_model), intent(in) :: this
    logical :: rst

    rst = this%mass_enabled
end function

subroutine dispatch_mass_matrix(this, x, y, m, args)
    class(mass_dispatch_model), intent(inout) :: this
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:,:) :: m
    class(*), intent(inout), optional :: args

    if (.not.this%mass_enabled) error stop "ODE path evaluated a mass matrix"
    this%mass_calls = this%mass_calls + 1
    m = 0.0d0
    m(1,1) = 2.0d0
    m(2,2) = 3.0d0
end subroutine

function check_mass_dispatch(integrator) result(rst)
    class(ode_integrator), intent(inout) :: integrator
    logical :: rst
    type(mass_dispatch_model) :: mdl
    integer :: mode
    real(real64), dimension(2) :: expected
    real(real64), allocatable, dimension(:,:) :: sol

    rst = .true.
    mdl%fcn => counted_mass_ode
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    do mode = 1, 4
        mdl%mass_enabled = mode == 2 .or. mode == 3
        mdl%mass_calls = 0
        call mdl%set_is_mass_matrix_dependent(mode == 3)
        call integrator%clear_buffer()
        call integrator%solve(mdl, [0.0d0, 0.1d0], [1.0d0, 1.0d0])
        sol = integrator%get_solution()
        if (mdl%mass_enabled) then
            expected = [exp(-0.1d0), exp(-0.2d0)]
            if (mode == 2) rst = rst .and. mdl%mass_calls == 1
            if (mode == 3) rst = rst .and. mdl%mass_calls > 1
        else
            expected = [exp(-0.2d0), exp(-0.6d0)]
            rst = rst .and. mdl%mass_calls == 0
        end if
        rst = rst .and. all(abs(sol(size(sol,1),2:3) - expected) < 1.0d-7)
    end do
end function

function test_solver_mass_dispatch() result(rst)
    logical :: rst
    logical, dimension(5) :: results
    type(bdf) :: bdf_integrator
    type(adams) :: adams_integrator
    type(rosenbrock) :: rosenbrock_integrator
    type(kennedy_carpenter_4) :: kc4_integrator
    type(kennedy_carpenter_5) :: kc5_integrator

    results(1) = check_mass_dispatch(bdf_integrator)
    results(2) = check_mass_dispatch(adams_integrator)
    results(3) = check_mass_dispatch(rosenbrock_integrator)
    results(4) = check_mass_dispatch(kc4_integrator)
    results(5) = check_mass_dispatch(kc5_integrator)
    rst = all(results)
    if (.not.rst) print *, "TEST FAILED: test_solver_mass_dispatch", results
end function

! ------------------------------------------------------------------------------
function test_analytical_jacobian_usage() result(rst)
    logical :: rst
    type(rosenbrock) :: integrator
    type(rosenbrock) :: mass_integrator
    type(ode_container) :: mdl
    type(jacobian_counter) :: counter
    real(real64), allocatable :: sol(:,:)

    mdl%fcn => counted_linear_ode
    mdl%jacobian => counted_linear_jacobian
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 1.0d0], counter)
    sol = integrator%get_solution()
    rst = counter%calls > 0 .and. &
        abs(sol(size(sol,1),2) - exp(-1.0d0)) < 1.0d-6 .and. &
        abs(sol(size(sol,1),3) - exp(-2.0d0)) < 1.0d-6

    counter%calls = 0
    counter%mass_calls = 0
    mdl%fcn => counted_mass_ode
    mdl%mass_matrix => counted_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call mass_integrator%set_absolute_tolerance(1.0d-9)
    call mass_integrator%set_relative_tolerance(1.0d-9)
    call mass_integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 1.0d0], counter)
    sol = mass_integrator%get_solution()
    rst = rst .and. counter%mass_calls == 1 .and. &
        abs(sol(size(sol,1),2) - exp(-1.0d0)) < 1.0d-6 .and. &
        abs(sol(size(sol,1),3) - exp(-2.0d0)) < 1.0d-6
end function

! ------------------------------------------------------------------------------
function test_stiff_vanderpol() result(rst)
    logical :: rst
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    real(real64) :: mu
    real(real64), allocatable :: sol(:,:)

    mu = 1.0d2
    mdl%fcn => vanderpol_args
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)
    call integrator%solve(mdl, [0.0d0, 2.0d0], [2.0d0, 0.0d0], mu)
    sol = integrator%get_solution()
    rst = size(sol, 1) > 2 .and. all(sol == sol)
end function

! ------------------------------------------------------------------------------
subroutine counted_linear_ode(x, y, dydx, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    dydx(1) = -y(1)
    dydx(2) = -2.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
subroutine counted_linear_jacobian(x, y, jac, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: jac(:,:)
    class(*), intent(inout), optional :: args

    select type (args)
    type is (jacobian_counter)
        args%calls = args%calls + 1
    end select
    jac = 0.0d0
    jac(1,1) = -1.0d0
    jac(2,2) = -2.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine counted_mass_ode(x, y, dydx, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    dydx(1) = -2.0d0 * y(1)
    dydx(2) = -6.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
subroutine counted_mass_matrix(x, y, mass, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: mass(:,:)
    class(*), intent(inout), optional :: args

    select type (args)
    type is (jacobian_counter)
        args%mass_calls = args%mass_calls + 1
    end select
    mass = 0.0d0
    mass(1,1) = 2.0d0
    mass(2,2) = 3.0d0
end subroutine

! ------------------------------------------------------------------------------
function test_implicit_rk_state_tolerances() result(rst)
    logical :: rst
    type(rosenbrock) :: rosenbrock_integrator
    type(kennedy_carpenter_4) :: kc_integrator
    type(ode_container) :: mdl
    real(real64) :: rosenbrock_error, kc_error
    real(real64), allocatable :: rosenbrock_sol(:,:), kc_sol(:,:), ans(:)

    call rosenbrock_integrator%set_absolute_tolerance([1.0d0, 2.0d0])
    call rosenbrock_integrator%set_relative_tolerance([1.0d-1, 2.0d-1])
    call kc_integrator%set_absolute_tolerance([1.0d0, 2.0d0])
    call kc_integrator%set_relative_tolerance([1.0d-1, 2.0d-1])
    rosenbrock_error = rosenbrock_integrator%compute_error_norm( &
        [1.0d1, 2.0d1], [1.0d1, 2.0d1], [2.0d0, 6.0d0])
    kc_error = kc_integrator%compute_error_norm([1.0d1, 2.0d1], &
        [1.0d1, 2.0d1], [2.0d0, 6.0d0])
    rst = abs(rosenbrock_error - 1.0d0) < epsilon(1.0d0) .and. &
        abs(kc_error - 1.0d0) < epsilon(1.0d0)

    mdl%fcn => test_2dof_1
    call rosenbrock_integrator%set_absolute_tolerance([1.0d-10, 1.0d-11])
    call rosenbrock_integrator%set_relative_tolerance([1.0d-11, 1.0d-10])
    call kc_integrator%set_absolute_tolerance([1.0d-10, 1.0d-11])
    call kc_integrator%set_relative_tolerance([1.0d-11, 1.0d-10])
    call rosenbrock_integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    call kc_integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    rosenbrock_sol = rosenbrock_integrator%get_solution()
    kc_sol = kc_integrator%get_solution()
    ans = test_2dof_solution_1(rosenbrock_sol(:,1))
    rst = rst .and. assert(ans, rosenbrock_sol(:,2), 1.0d-5)
    ans = test_2dof_solution_1(kc_sol(:,1))
    rst = rst .and. assert(ans, kc_sol(:,2), 1.0d-5)
end function

! ------------------------------------------------------------------------------
function test_kennedy_carpenter_4() result(rst)
    logical :: rst
    type(kennedy_carpenter_4) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    ans = test_1dof_solution_1(sol(:,1))
    rst = assert(ans, sol(:,2), 1.0d-5)
end function

! ------------------------------------------------------------------------------
function test_kennedy_carpenter_5() result(rst)
    logical :: rst
    type(kennedy_carpenter_5) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    ans = test_1dof_solution_1(sol(:,1))
    rst = assert(ans, sol(:,2), 1.0d-6)
end function

! ------------------------------------------------------------------------------
function test_kennedy_carpenter_newton_settings() result(rst)
    logical :: rst
    type(kennedy_carpenter_4) :: kc4
    type(kennedy_carpenter_5) :: kc5
    type(ode_container) :: mdl
    type(ode_results) :: res
    real(real64), allocatable :: sol(:,:), ans(:)
    integer(int32) :: nsteps

    rst = .true.

    ! Defaults
    if (kc4%get_newton_step_limit() /= 12 .or. &
        kc5%get_newton_step_limit() /= 12) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter default Newton step limit"
    end if
    if (abs(kc4%get_newton_tolerance() - 1.0d-8) > 1.0d-20 .or. &
        abs(kc5%get_newton_tolerance() - 1.0d-8) > 1.0d-20) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter default Newton tolerance"
    end if

    ! Round trip
    call kc4%set_newton_step_limit(20)
    call kc4%set_newton_tolerance(1.0d-11)
    call kc5%set_newton_step_limit(5)
    call kc5%set_newton_tolerance(1.0d-6)
    if (kc4%get_newton_step_limit() /= 20 .or. &
        kc5%get_newton_step_limit() /= 5) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter Newton step limit set/get"
    end if
    if (abs(kc4%get_newton_tolerance() - 1.0d-11) > 1.0d-24 .or. &
        abs(kc5%get_newton_tolerance() - 1.0d-6) > 1.0d-20) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter Newton tolerance set/get"
    end if

    ! Solve with the modified settings and verify accuracy and that the
    ! Jacobian is reused across steps
    mdl%fcn => test_1dof_1
    call kc4%set_absolute_tolerance(1.0d-10)
    call kc4%set_relative_tolerance(1.0d-10)
    call kc4%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    res = kc4%get_ode_results()
    sol = kc4%get_solution()
    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), 1.0d-5)) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter solution w/ Newton settings"
    end if
    nsteps = res%statistics%successful_steps + res%statistics%failed_steps
    if (res%statistics%jacobian_evaluations >= nsteps) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter Jacobian reuse"
    end if

    call kc5%set_absolute_tolerance(1.0d-8)
    call kc5%set_relative_tolerance(1.0d-8)
    call kc5%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = kc5%get_solution()
    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), 1.0d-5)) then
        rst = .false.
        print "(A)", "TEST FAILED: Kennedy-Carpenter 5 solution w/ Newton settings"
    end if
end function

! ------------------------------------------------------------------------------
function test_kennedy_carpenter_mass_matrix() result(rst)
    logical :: rst
    type(kennedy_carpenter_4) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)
    mdl%fcn => rosenbrock_mass_ode
    mdl%mass_matrix => rosenbrock_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 1.0d0])
    sol = integrator%get_solution()
    rst = abs(sol(size(sol,1),2) - exp(-1.0d0)) < 1.0d-7 .and. &
        abs(sol(size(sol,1),3) - exp(-2.0d0)) < 1.0d-7
end function

! ------------------------------------------------------------------------------
function test_kennedy_carpenter_singular_mass_matrix() result(rst)
    logical :: rst
    type(kennedy_carpenter_4) :: integrator4
    type(kennedy_carpenter_5) :: integrator5
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable :: sol4(:,:), sol5(:,:)
    real(real64) :: length

    length = 1.5d0
    args%length = length
    args%mass = 2.0d0
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator4%set_absolute_tolerance(1.0d-9)
    call integrator4%set_relative_tolerance(1.0d-9)
    call integrator5%set_absolute_tolerance(1.0d-9)
    call integrator5%set_relative_tolerance(1.0d-9)

    call integrator4%solve(mdl, [0.0d0, 0.1d0], &
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], args)
    call integrator5%solve(mdl, [0.0d0, 0.1d0], &
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], args)
    sol4 = integrator4%get_solution()
    sol5 = integrator5%get_solution()

    rst = abs(sol4(size(sol4,1),2)**2 + sol4(size(sol4,1),4)**2 - &
        length**2) < 1.0d-7 .and. &
        abs(sol5(size(sol5,1),2)**2 + sol5(size(sol5,1),4)**2 - &
        length**2) < 1.0d-7
end function

! ------------------------------------------------------------------------------
function test_rosenbrock_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_1dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_rosenbrock_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_rosenbrock_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    integer :: i
    real(real64) :: x(npts)
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1

    ! Define the values where to compute the solution
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Compute the solution
    call integrator%solve(mdl, x, [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_3 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_rosenbrock_mass_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5
    integer, parameter :: npts = 1000
    real(real64), parameter :: tmax = 1.0d1

    ! Local Variables
    integer :: i
    real(real64) :: dt, t(npts)
    type(rosenbrock) :: integrator
    type(ode_container) :: mass_mdl, ref_mdl
    real(real64), allocatable, dimension(:,:) :: sol, refsol

    ! Initialization
    rst = .true.
    dt = tmax / (npts - 1.0d0)
    t = (/ (i * dt, i = 0, npts - 1) /)
    mass_mdl%fcn => rosenbrock_mass_ode
    mass_mdl%mass_matrix => rosenbrock_mass_matrix
    call mass_mdl%set_is_mass_matrix_dependent(.false.)
    ref_mdl%fcn => rosenbrock_reference_ode

    ! Solve the mass-matrix form and the equivalent standard ODE form
    call integrator%solve(mass_mdl, t, [1.0d0, 1.0d0])
    sol = integrator%get_solution()

    call integrator%clear_buffer()
    call integrator%solve(ref_mdl, t, [1.0d0, 1.0d0])
    refsol = integrator%get_solution()

    ! Test
    if (.not.assert(sol, refsol, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_mass_matrix -1"
    end if

    call integrator%clear_buffer()
    call integrator%solve(mass_mdl, t, [1.0d0, 1.0d0])
    sol = integrator%get_solution()
    if (.not.assert(sol, refsol, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_mass_matrix -2"
    end if
end function

! ------------------------------------------------------------------------------
subroutine rosenbrock_mass_matrix(x, y, m, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:,:) :: m
    class(*), intent(inout), optional :: args

    m(1,1) = 2.0d0
    m(1,2) = 0.0d0
    m(2,1) = 0.0d0
    m(2,2) = 3.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine rosenbrock_mass_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    ! M * y' = g(y), with M = diag(2,3) and g(y) = [-2*y1, -6*y2]
    ! This is equivalent to y' = [-y1, -2*y2].
    dydx(1) = -2.0d0 * y(1)
    dydx(2) = -6.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
subroutine rosenbrock_reference_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    dydx(1) = -1.0d0 * y(1)
    dydx(2) = -2.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
function test_rosenbrock_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3
    integer, parameter :: npts = 1000
    real(real64), parameter :: tmax = 5.0d1

    ! Local Variables
    integer :: i
    real(real64) :: mu, dt, t(npts)
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl, ref
    real(real64), allocatable, dimension(:,:) :: sol, refsol

    ! Initialization
    rst = .true.
    mdl%fcn => vanderpol_args
    ref%fcn => vanderpol
    mu = 5.0d0
    dt = tmax / (npts - 1.0d0)
    t = (/ (i * dt, i = 0, npts - 1) /)

    ! Perform the integration with user-defined arguments
    call integrator%solve(mdl, t, [2.0d0, 0.0d0], args = mu)
    sol = integrator%get_solution()

    ! Perform the integration without additional arguments
    call integrator%clear_buffer()
    call integrator%solve(ref, t, [2.0d0, 0.0d0])
    refsol = integrator%get_solution()

    ! Test
    if (.not.assert(sol, refsol, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_with_args -1"
        print 100, "Solution Size: ", size(sol, 1), "-", size(sol, 2)
        print 100, "Reference Size: ", size(refsol, 1), "-", size(refsol, 2)
        print 101, "Solution - Reference Norm 1: ", norm2(sol(:,2) - refsol(:,2))
        print 101, "Solution - Reference Norm 2: ", norm2(sol(:,3) - refsol(:,3))
        print 101, "Max Delta 1: ", maxval(abs(sol(:,2) - refsol(:,2)))
        print 101, "Max Delta 2: ", maxval(abs(sol(:,3) - refsol(:,3)))
    end if

    ! Formatting
100 format(A, I0, A, I0)
101 format(A, G12.3)
end function

! ------------------------------------------------------------------------------
end module