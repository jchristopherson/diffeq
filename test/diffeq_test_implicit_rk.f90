module diffeq_test_implicit_rk
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
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
function test_rosenbrock_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: mu
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl, ref
    real(real64), allocatable, dimension(:,:) :: sol, refsol

    ! Initialization
    rst = .true.
    mdl%fcn => vanderpol_args
    ref%fcn => vanderpol
    mu = 5.0d0

    ! Perform the integration with user-defined arguments
    call integrator%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0], args = mu)
    sol = integrator%get_solution()

    ! Perform the integration without additional arguments
    call integrator%clear_buffer()
    call integrator%solve(ref, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
    refsol = integrator%get_solution()

    ! Test
    if (.not.assert(sol, refsol, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rosenbrock_with_args -1"
        print 100, "Solution Size: ", size(sol, 1), "-", size(sol, 2)
        print 100, "Reference Size: ", size(ref, 1), "-", size(ref, 2)
    end if

    ! Formatting
100 format(A, I0, A, I0)
end function

! ------------------------------------------------------------------------------
end module