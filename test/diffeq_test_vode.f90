module diffeq_test_vode
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
! ------------------------------------------------------------------------------
function test_adams_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tmax = 1.0d0
    real(real64), parameter :: ic(2) = [1.0d0, 0.5d0]
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: ans, t
    real(real64), allocatable, dimension(:,:) :: sol

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, tmax], ic)
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_1 -1"
    end if

    ! Test 2
    allocate(t(npts))
    dt = tmax / (npts - 1.0d0)
    t = (/ (i * dt, i = 0, npts - 1) /)
    
    call integrator%clear_buffer()
    call integrator%solve(mdl, t, ic)
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_bdf_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tmax = 1.0d0
    real(real64), parameter :: ic(2) = [1.0d0, 0.5d0]
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(bdf) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: ans, t
    real(real64), allocatable, dimension(:,:) :: sol

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, tmax], ic)
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_1 -1"
    end if

    ! Test 2
    allocate(t(npts))
    dt = tmax / (npts - 1.0d0)
    t = (/ (i * dt, i = 0, npts - 1) /)
    
    call integrator%clear_buffer()
    call integrator%solve(mdl, t, ic)
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_1 -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-2
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tmax = 5.0d1

    ! Local Variables
    integer(int32) :: i
    real(real64) :: mu, dt, t(npts)
    type(adams) :: integrator
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
        print "(A)", "TEST FAILED: test_adams_with_args -1"
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
function test_bdf_with_args() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-2
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: tmax = 5.0d1

    ! Local Variables
    integer(int32) :: i
    real(real64) :: mu, dt, t(npts)
    type(adams) :: integrator
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
        print "(A)", "TEST FAILED: test_bdf_with_args -1"
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