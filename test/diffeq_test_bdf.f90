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

module diffeq_test_bdf
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
function test_bdf_nonlinear_reaction() result(rst)
    logical :: rst
    type(bdf) :: integrator
    type(ode_container) :: model
    integer, parameter :: point_count = 301
    real(real64), dimension(point_count) :: grid, expected
    real(real64), allocatable, dimension(:,:) :: solution
    integer :: index

    rst = .true.
    model%fcn => nonlinear_reaction_ode
    do index = 1, point_count
        grid(index) = (real(index - 1, real64) / real(point_count - 1, real64))**2
    end do
    call integrator%set_absolute_tolerance([1.0d-10, 1.0d-12, 1.0d-10])
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%set_step_limit(100000)
    call integrator%solve(model, grid, [1.0d-3, 9.99d-1, 1.0d0])
    solution = integrator%get_solution()
    if (.not.all(shape(solution) == [point_count, 4])) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_nonlinear_reaction shape"
        return
    end if
    expected = 1.0d0 / (1.0d0 + 999.0d0 * exp(-1000.0d0 * grid))
    if (.not.all(ieee_is_finite(solution)) .or. &
        .not.assert(grid, solution(:,1), 1.0d-12) .or. &
        .not.assert(expected, solution(:,2), 5.0d-7) .or. &
        .not.assert(1.0d0 - expected, solution(:,3), 5.0d-7) .or. &
        .not.assert(exp(-grid), solution(:,4), 5.0d-7) .or. &
        maxval(abs(solution(:,2) + solution(:,3) - 1.0d0)) > 1.0d-8 .or. &
        any(solution(:,2:) < -1.0d-8)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_nonlinear_reaction accuracy or invariants"
    end if
end function

subroutine nonlinear_reaction_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    dydx(1) = 1000.0d0 * y(1) * y(2)
    dydx(2) = -dydx(1)
    dydx(3) = -y(3)
end subroutine

function test_bdf_robustness() result(rst)
    logical :: rst
    type(bdf) :: integrator

    rst = check_multistep_robustness(integrator, "BDF")
end function

function check_multistep_robustness(integrator, name) result(rst)
    class(ode_integrator), intent(inout) :: integrator
    character(len=*), intent(in) :: name
    logical :: rst
    type(ode_container) :: model
    integer, parameter :: point_count = 2501
    real(real64), dimension(point_count) :: grid, expected_position, expected_velocity
    real(real64), dimension(2) :: initial
    real(real64), allocatable, dimension(:,:) :: solution
    real(real64) :: start, finish
    integer :: index, direction, count

    rst = .true.
    model%fcn => test_2dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    call integrator%set_step_limit(100000)
    do direction = 1, 2
        start = 0.0d0
        finish = 20.0d0
        if (direction == 2) then
            start = 20.0d0
            finish = 0.0d0
        end if
        do index = 1, point_count
            grid(index) = start + (finish - start) * &
                (real(index - 1, real64) / real(point_count - 1, real64))**2
        end do
        call test_2dof_state_solution_1(start, initial)
        call integrator%clear_buffer()
        call integrator%set_allow_overshoot(direction == 1)
        call integrator%solve(model, grid, initial)
        solution = integrator%get_solution()
        if (.not.all(shape(solution) == [point_count, 3])) then
            rst = .false.
            print "(A)", "TEST FAILED: " // name // " long-run shape"
            return
        end if
        expected_position = test_2dof_solution_1(grid)
        expected_velocity = 0.5d0 * cos(20.0d0 * grid) - 20.0d0 * sin(20.0d0 * grid)
        if (.not.all(ieee_is_finite(solution)) .or. &
            .not.assert(grid, solution(:,1), 1.0d-12) .or. &
            .not.assert(expected_position, solution(:,2), 2.0d-6) .or. &
            .not.assert(expected_velocity, solution(:,3), 4.0d-5)) then
            rst = .false.
            print "(A,I0)", "TEST FAILED: " // name // " long-run direction ", direction
        end if
    end do

    call integrator%clear_buffer()
    model%fcn => robustness_decay_ode
    call integrator%set_absolute_tolerance([1.0d-2, 1.0d-18])
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%solve(model, [0.0d0, 0.1d0, 0.3d0, 1.0d0, 2.0d0], [1.0d8, 1.0d-8])
    solution = integrator%get_solution()
    if (.not.all(shape(solution) == [5, 3])) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " mixed-scale shape"
        return
    end if
    if (.not.all(ieee_is_finite(solution)) .or. &
        .not.assert(exp(-solution(:,1)), solution(:,2) / 1.0d8, 2.0d-7) .or. &
        .not.assert(exp(-solution(:,1)), solution(:,3) / 1.0d-8, 2.0d-7)) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " mixed-scale accuracy"
    end if

    call integrator%clear_buffer()
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%solve(model, [0.0d0, 0.1d0, 1.0d0], [0.0d0, 0.0d0, 0.0d0])
    solution = integrator%get_solution()
    if (.not.all(shape(solution) == [3, 4])) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " equilibrium shape"
        return
    end if
    if (.not.all(ieee_is_finite(solution)) .or. any(solution(:,2:) /= 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " equilibrium"
    end if

    call integrator%clear_buffer()
    call integrator%set_maximum_step_size(5.0d-4)
    call integrator%solve(model, [0.0d0, 1.0d0], [1.0d0])
    solution = integrator%get_solution()
    count = size(solution, 1)
    if (count <= 1000 .or. size(solution, 2) /= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " accepted-step buffer growth"
        return
    end if
    if (.not.all(ieee_is_finite(solution)) .or. &
        any(solution(2:,1) <= solution(:count - 1,1)) .or. &
        maxval(solution(2:,1) - solution(:count - 1,1)) > 5.0d-4 * (1.0d0 + 1.0d-10) .or. &
        abs(solution(1,1)) > 1.0d-12 .or. abs(solution(count,1) - 1.0d0) > 1.0d-12 .or. &
        .not.assert(exp(-solution(:,1)), solution(:,2), 2.0d-7)) then
        rst = .false.
        print "(A)", "TEST FAILED: " // name // " accepted-step accuracy or limits"
    end if
end function

subroutine robustness_decay_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args
    dydx = -y
end subroutine

function test_bdf_stiff_relaxation() result(rst)
    logical :: rst
    type(bdf) :: integrator
    type(ode_container) :: model
    integer, parameter :: point_count = 251
    real(real64), dimension(point_count) :: grid, expected
    real(real64), dimension(3) :: rates
    real(real64), allocatable, dimension(:,:) :: solution
    real(real64) :: rate
    integer :: index, trial

    rst = .true.
    model%fcn => stiff_relaxation_ode
    rates = [1.0d1, 1.0d3, 1.0d5]
    do index = 1, point_count
        grid(index) = (real(index - 1, real64) / real(point_count - 1, real64))**2
    end do
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%set_step_limit(100000)
    do trial = 1, size(rates)
        rate = rates(trial)
        call integrator%clear_buffer()
        call integrator%solve(model, grid, [2.0d0], rate)
        solution = integrator%get_solution()
        if (.not.all(shape(solution) == [point_count, 2])) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bdf_stiff_relaxation shape"
            return
        end if
        expected = cos(grid) + exp(-rate * grid)
        if (.not.all(ieee_is_finite(solution)) .or. &
            .not.assert(grid, solution(:,1), 1.0d-12) .or. &
            .not.assert(expected, solution(:,2), 5.0d-7)) then
            rst = .false.
            print "(A,ES12.4)", "TEST FAILED: test_bdf_stiff_relaxation rate ", rate
        end if
    end do
end function

subroutine stiff_relaxation_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    select type (args)
    type is (real(real64))
        dydx(1) = -args * (y(1) - cos(x)) - sin(x)
    class default
        error stop "Missing relaxation rate"
    end select
end subroutine

! ------------------------------------------------------------------------------
function test_bdf_reverse() result(rst)
    logical :: rst
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)

    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    call integrator%solve(mdl, [1.0d0, 0.0d0], &
        [test_1dof_solution_1(1.0d0)])
    sol = integrator%get_solution()
    rst = abs(sol(size(sol,1),1)) < 1.0d-12 .and. &
        abs(sol(size(sol,1),2) - 2.0d0) < 1.0d-7
end function

! ------------------------------------------------------------------------------
function test_bdf_state_tolerances() result(rst)
    logical :: rst
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64) :: err
    real(real64), allocatable :: sol(:,:), ans(:)

    call integrator%set_absolute_tolerance([1.0d0, 2.0d0])
    call integrator%set_relative_tolerance([1.0d-1, 2.0d-1])
    err = integrator%compute_error_norm([1.0d1, 2.0d1], &
        [1.0d1, 2.0d1], [2.0d0, 6.0d0])
    rst = abs(err - 1.0d0) < epsilon(1.0d0)

    mdl%fcn => test_2dof_1
    call integrator%set_absolute_tolerance([1.0d-10, 1.0d-11])
    call integrator%set_relative_tolerance([1.0d-11, 1.0d-10])
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    sol = integrator%get_solution()
    ans = test_2dof_solution_1(sol(:,1))
    rst = rst .and. assert(ans, sol(:,2), 1.0d-5)
end function

! ------------------------------------------------------------------------------
function test_bdf_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_1dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_dense() result(rst)
    ! Tests the dense output capability by requesting the solution at a set of
    ! points that bear no relationship to the steps the solver actually takes.

    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-3
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    integer :: i
    real(real64) :: x(npts)
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Compute the solution
    call integrator%solve(mdl, x, [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! The solution must be returned at precisely the requested points
    if (size(sol, 1) /= npts) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_dense -1"
        return
    end if
    if (.not.assert(x, sol(:,1), 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_dense -2"
    end if

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_dense -3"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_all_steps() result(rst)
    ! Tests that supplying only the end points of the integration range causes
    ! every successful step to be returned.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    integer :: i, n
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    n = size(sol, 1)

    ! More than just the end points must be reported
    if (n <= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_all_steps -1"
        return
    end if

    ! The steps must advance monotonically and span the requested range
    do i = 2, n
        if (sol(i,1) <= sol(i-1,1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bdf_all_steps -2"
            return
        end if
    end do
    if (abs(sol(1,1)) > 1.0d-12 .or. abs(sol(n,1) - 1.0d0) > 1.0d-12) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_all_steps -3"
    end if

    ! Every reported step must be accurate
    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_all_steps -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_mass_matrix() result(rst)
    ! Tests a system carrying a nonsingular mass matrix against the equivalent
    ! system written in explicit form.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)
    integer :: n

    ! Initialization
    rst = .true.
    mdl%fcn => bdf_mass_ode
    mdl%mass_matrix => bdf_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Compute the solution
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 1.0d0])
    sol = integrator%get_solution()
    n = size(sol, 1)

    ! M * y' = g(y) with M = diag(2,3) reduces to y' = [-y1, -2*y2]
    if (abs(sol(n,2) - exp(-1.0d0)) > tol .or. &
        abs(sol(n,3) - exp(-2.0d0)) > tol) &
    then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_mass_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
subroutine bdf_mass_matrix(x, y, m, args)
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
subroutine bdf_mass_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    ! M * y' = g(y), with M = diag(2,3) and g(y) = [-2*y1, -6*y2]
    dydx(1) = -2.0d0 * y(1)
    dydx(2) = -6.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
function test_bdf_singular_mass_matrix() result(rst)
    ! Tests an index-1 DAE.  The Cartesian pendulum carries a mass matrix whose
    ! final row is identically zero, so the corresponding equation is an
    ! algebraic constraint rather than a differential one.  A correct solution
    ! keeps the bob on the circle of radius L.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable :: sol(:,:)
    integer :: i, n
    real(real64) :: r

    ! Initialization
    rst = .true.
    args%length = length
    args%mass = 2.0d0
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)

    ! Compute the solution
    call integrator%solve(mdl, [0.0d0, 1.0d0], &
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], args)
    sol = integrator%get_solution()
    n = size(sol, 1)

    ! The constraint must hold at every reported point
    do i = 1, n
        r = sol(i,2)**2 + sol(i,4)**2 - length**2
        if (abs(r) > tol) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bdf_singular_mass_matrix -1"
            return
        end if
    end do

    ! The pendulum must actually have moved
    if (abs(sol(n,4)) < 1.0d-3) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_singular_mass_matrix -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_singular_mass_matrix_dense() result(rst)
    ! Tests dense output for an index-1 DAE.  Interpolated points must satisfy
    ! the algebraic constraint just as the accepted steps do.

    ! Arguments
    logical :: rst

    ! Parameters
    integer, parameter :: npts = 200
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable :: sol(:,:)
    integer :: i
    real(real64) :: dt, t(npts), r

    ! Initialization
    rst = .true.
    args%length = length
    args%mass = 2.0d0
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    dt = 1.0d0 / (npts - 1.0d0)
    t = (/ (i * dt, i = 0, npts - 1) /)

    ! Compute the solution
    call integrator%solve(mdl, t, &
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], args)
    sol = integrator%get_solution()

    ! The solution must be returned at the requested points
    if (size(sol, 1) /= npts) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_singular_mass_matrix_dense -1"
        return
    end if

    ! The constraint must hold at every interpolated point
    do i = 1, npts
        r = sol(i,2)**2 + sol(i,4)**2 - length**2
        if (abs(r) > tol) then
            rst = .false.
            print "(A)", "TEST FAILED: test_bdf_singular_mass_matrix_dense -2"
            return
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_bdf_with_args() result(rst)
    ! Tests that user-defined arguments reach the model.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl, refmdl
    real(real64) :: mu
    real(real64), allocatable :: sol(:,:), ref(:,:)
    integer :: i
    real(real64) :: x(100)

    ! Initialization
    rst = .true.
    mu = 5.0d0
    x = (/ (i * 1.0d-1, i = 0, 99) /)
    mdl%fcn => vanderpol_args
    refmdl%fcn => vanderpol
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)

    ! Solve with the parameter supplied through args
    call integrator%solve(mdl, x, [2.0d0, 0.0d0], mu)
    sol = integrator%get_solution()

    ! Solve the equivalent model carrying the parameter internally
    call integrator%clear_buffer()
    call integrator%solve(refmdl, x, [2.0d0, 0.0d0])
    ref = integrator%get_solution()

    ! Test
    if (.not.assert(sol, ref, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_with_args -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_order_range() result(rst)
    ! Tests that the order remains within the range the method supports and
    ! that the order limit may be constrained by the caller.

    ! Arguments
    logical :: rst

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1

    ! The method must report its documented order limit
    if (integrator%get_maximum_supported_order() /= 5) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_order_range -1"
    end if

    ! A request beyond that limit must be clamped
    call integrator%set_maximum_order(9)
    if (integrator%get_maximum_order() /= 5) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_order_range -2"
    end if

    ! Restricting the order must still produce a usable solution
    call integrator%set_maximum_order(2)
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()

    if (integrator%get_order() > 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_order_range -3"
    end if

    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), 1.0d-4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_order_range -4"
    end if
end function

! ------------------------------------------------------------------------------
end module
