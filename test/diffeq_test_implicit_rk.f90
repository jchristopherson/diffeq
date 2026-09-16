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

contains
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
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
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
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tmax = 1.0d1

    ! Local Variables
    integer(int32) :: i
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
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tmax = 5.0d1

    ! Local Variables
    integer(int32) :: i
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