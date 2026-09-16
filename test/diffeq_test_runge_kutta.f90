module diffeq_test_runge_kutta
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
! ------------------------------------------------------------------------------
function test_reverse_and_overshoot() result(rst)
    logical :: rst
    type(runge_kutta_45) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)

    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    call integrator%solve(mdl, [1.0d0, 0.0d0], &
        [test_1dof_solution_1(1.0d0)])
    sol = integrator%get_solution()
    rst = abs(sol(size(sol,1),1)) < 1.0d-12 .and. &
        abs(sol(size(sol,1),2) - 2.0d0) < 1.0d-8

    call integrator%clear_buffer()
    call integrator%set_absolute_tolerance(1.0d0)
    call integrator%set_relative_tolerance(1.0d0)
    call integrator%set_maximum_step_size(2.0d-1)
    call integrator%set_allow_overshoot(.true.)
    call integrator%solve(mdl, [0.0d0, 3.1d-1], [2.0d0])
    sol = integrator%get_solution()
    rst = rst .and. abs(sol(size(sol,1),1) - 3.1d-1) < 1.0d-12

    call integrator%clear_buffer()
    call integrator%set_allow_overshoot(.false.)
    call integrator%solve(mdl, [0.0d0, 1.0d-1, 2.0d-1, 3.1d-1], [2.0d0])
    sol = integrator%get_solution()
    rst = rst .and. abs(sol(size(sol,1),1) - 3.1d-1) < 1.0d-12 .and. &
        size(sol,1) == 4
end function

! ------------------------------------------------------------------------------
function test_state_variable_tolerances() result(rst)
    logical :: rst
    type(tsitouras_54) :: integrator
    type(runge_kutta_853) :: integrator853
    type(ode_container) :: mdl
    real(real64) :: err
    real(real64), allocatable :: sol(:,:), ans(:)

    call integrator%set_absolute_tolerance(1.0d-3)
    call integrator%set_relative_tolerance(2.0d-2)
    rst = integrator%get_absolute_tolerance() == 1.0d-3 .and. &
        integrator%get_relative_tolerance() == 2.0d-2

    call integrator%set_absolute_tolerance([1.0d0, 2.0d0])
    call integrator%set_relative_tolerance([1.0d-1, 2.0d-1])
    err = integrator%compute_error_norm([1.0d1, 2.0d1], &
        [1.0d1, 2.0d1], [2.0d0, 6.0d0])
    rst = rst .and. abs(err - 1.0d0) < epsilon(1.0d0) .and. &
        integrator%get_absolute_tolerance(1) == 1.0d0 .and. &
        integrator%get_relative_tolerance(2) == 2.0d-1

    call integrator%set_absolute_tolerance(2.0d0)
    call integrator%set_relative_tolerance(0.0d0)
    err = integrator%compute_error_norm([1.0d1, 2.0d1], &
        [1.0d1, 2.0d1], [2.0d0, 2.0d0])
    rst = rst .and. abs(err - 1.0d0) < epsilon(1.0d0)

    mdl%fcn => test_2dof_1
    call integrator853%set_absolute_tolerance([1.0d-7, 1.0d-6])
    call integrator853%set_relative_tolerance([1.0d-6, 1.0d-7])
    call integrator853%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    sol = integrator853%get_solution()
    ans = test_2dof_solution_1(sol(:,1))
    rst = rst .and. assert(ans, sol(:,2), 1.0d-5)
end function

! ------------------------------------------------------------------------------
function test_step_size_limits() result(rst)
    logical :: rst
    type(runge_kutta_45) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)
    real(real64) :: max_step

    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d0)
    call integrator%set_relative_tolerance(1.0d0)
    call integrator%set_maximum_step_size(5.0d-2)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    max_step = maxval(abs(sol(2:,1) - sol(:size(sol,1)-1,1)))
    rst = max_step <= 5.0d-2 * (1.0d0 + 10.0d0 * epsilon(1.0d0))
end function

! ------------------------------------------------------------------------------
function test_tsitouras_54() result(rst)
    logical :: rst
    type(tsitouras_54) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-9)
    call integrator%set_relative_tolerance(1.0d-9)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    ans = test_1dof_solution_1(sol(:,1))
    rst = integrator%get_order() == 5 .and. integrator%get_is_fsal() .and. &
        integrator%get_stage_count() == 7 .and. assert(ans, sol(:,2), 1.0d-5)
end function

! ------------------------------------------------------------------------------
function test_tsitouras_54_dense() result(rst)
    logical :: rst
    type(tsitouras_54) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)
    real(real64) :: state(2), x(101)
    integer :: i
    mdl%fcn => test_2dof_1
    x = [(real(i, real64) / 100.0d0, i = 0, 100)]
    call integrator%solve(mdl, x, [1.0d0, 0.5d0])
    sol = integrator%get_solution()
    rst = size(sol, 1) == size(x) .and. &
        assert(x, sol(:,1), 1.0d-12)
    allocate(ans(size(sol,1)))
    ans = test_2dof_solution_1(sol(:,1))
    rst = rst .and. assert(ans, sol(:,2), 1.0d-4)
    do i = 1, size(sol, 1)
        call test_2dof_state_solution_1(sol(i,1), state)
        if (abs(sol(i,3) - state(2)) > 1.0d-3) rst = .false.
    end do
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_45_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(runge_kutta_45) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_45_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(runge_kutta_45) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_45_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    type(runge_kutta_45) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_3 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_23_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    type(runge_kutta_23) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_23_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_23_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 2.0d-2

    ! Local Variables
    type(runge_kutta_23) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_23_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_23_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-2

    ! Local Variables
    type(runge_kutta_23) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_23_3 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_853_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(runge_kutta_853) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_853_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(runge_kutta_853) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_853_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(runge_kutta_853) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_853_3 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: mu
    type(runge_kutta_23) :: rk23
    type(runge_kutta_45) :: rk45
    type(runge_kutta_853) :: rk853
    type(ode_container) :: mdl, ref
    real(real64), allocatable, dimension(:,:) :: sol23, sol45, sol853, &
        ref23, ref45, ref853

    ! Initialization
    rst = .true.
    mdl%fcn => vanderpol_args
    ref%fcn => vanderpol
    mu = 5.0d0

    ! Perform the integration with user-defined arguments
    call rk23%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0], args = mu)
    call rk45%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0], args = mu)
    call rk853%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0], args = mu)

    sol23 = rk23%get_solution()
    sol45 = rk45%get_solution()
    sol853 = rk853%get_solution()

    ! Perform the integration without additional arguments
    call rk23%clear_buffer()
    call rk45%clear_buffer()
    call rk853%clear_buffer()
    call rk23%solve(ref, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
    call rk45%solve(ref, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
    call rk853%solve(ref, [0.0d0, 5.0d1], [2.0d0, 0.0d0])

    ref23 = rk23%get_solution()
    ref45 = rk45%get_solution()
    ref853 = rk853%get_solution()

    ! Test
    if (.not.assert(sol23, ref23, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_with_args -1"
    end if

    if (.not.assert(sol45, ref45, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_with_args -2"
    end if

    if (.not.assert(sol853, ref853, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_with_args -3"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_dense_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6
    real(real64), parameter :: maxX = 5.0d1
    integer(int32), parameter :: npts = 1000

    ! Local Variables
    integer(int32) :: i
    real(real64) :: mu, dx
    real(real64), allocatable, dimension(:) :: x
    type(runge_kutta_23) :: rk23
    type(runge_kutta_45) :: rk45
    type(runge_kutta_853) :: rk853
    type(ode_container) :: mdl, ref
    real(real64), allocatable, dimension(:,:) :: sol23, sol45, sol853, &
        ref23, ref45, ref853

    ! Initialization
    rst = .true.
    mdl%fcn => vanderpol_args
    ref%fcn => vanderpol
    mu = 5.0d0
    allocate(x(npts))
    dx = maxX / (npts - 1.0d0)
    x = (/ (i * dx, i = 0, npts - 1) /)

    ! Perform the integration with user-defined arguments
    call rk23%solve(mdl, x, [2.0d0, 0.0d0], args = mu)
    call rk45%solve(mdl, x, [2.0d0, 0.0d0], args = mu)
    call rk853%solve(mdl, x, [2.0d0, 0.0d0], args = mu)

    sol23 = rk23%get_solution()
    sol45 = rk45%get_solution()
    sol853 = rk853%get_solution()

    ! Perform the integration without additional arguments
    call rk23%clear_buffer()
    call rk45%clear_buffer()
    call rk853%clear_buffer()
    call rk23%solve(ref, x, [2.0d0, 0.0d0])
    call rk45%solve(ref, x, [2.0d0, 0.0d0])
    call rk853%solve(ref, x, [2.0d0, 0.0d0])

    ref23 = rk23%get_solution()
    ref45 = rk45%get_solution()
    ref853 = rk853%get_solution()

    ! Test
    if (.not.assert(sol23, ref23, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_dense_with_args -1"
    end if

    if (.not.assert(sol45, ref45, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_dense_with_args -2"
    end if

    if (.not.assert(sol853, ref853, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_runge_kutta_dense_with_args -3"
    end if
end function

! ------------------------------------------------------------------------------
end module