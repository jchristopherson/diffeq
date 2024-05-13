module diffeq_test_runge_kutta
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_3 -1"
    end if
end function

! ------------------------------------------------------------------------------
end module