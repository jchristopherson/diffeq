module diffeq_sdirk4_tests
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_sdirk4_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    type(sdirk4_integrator) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.

    ! Define the model
    mdl%fcn => test_1dof_1

    ! Perform the integration
    sol = integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])

    ! Compute the actual solution
    ans = test_1dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_sdirk4_1 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_sdirk4_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 2.0d-3

    ! Local Variables
    type(sdirk4_integrator) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.

    ! Define the model
    mdl%fcn => test_2dof_1

    ! Perform the integration
    sol = integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_sdirk4_2 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_sdirk4_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(sdirk4_integrator) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts), ans(npts), sol(npts,3)

    ! Initialization
    rst = .true.

    ! Define the model
    mdl%fcn => test_2dof_1

    ! Define the values of x where the solution should be computed
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Compute the actual solution
    ans = test_2dof_solution_1(x)

    ! Compute the solution
    sol = integrator%solve(mdl, x, [1.0d0, 0.5d0])

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_sdirk4_3 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
end module