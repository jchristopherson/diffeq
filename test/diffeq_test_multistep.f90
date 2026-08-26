module diffeq_test_multistep
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
function test_multistep_method_contract() result(rst)
    logical :: rst
    type(adams) :: adams_solver
    type(bdf) :: bdf_solver
    type(ode_container) :: model
    real(real64), allocatable :: solution(:,:)
    rst = adams_solver%get_method() == DIFFEQ_ADAMS_METHOD .and. &
        bdf_solver%get_method() == DIFFEQ_BDF_METHOD .and. &
        adams_solver%get_order() == 12 .and. bdf_solver%get_order() == 5
    model%fcn => test_1dof_1
    call adams_solver%solve(model, [0.0d0, 1.0d0], [2.0d0])
    solution = adams_solver%get_solution()
    call bdf_solver%solve(model, [0.0d0, 1.0d0], [2.0d0])
    solution = bdf_solver%get_solution()
end function

function test_adams_high_order() result(rst)
    logical :: rst
    type(adams) :: solver
    type(ode_container) :: model
    real(real64), allocatable :: solution(:,:)
    model%fcn => exponential_model
    call solver%solve(model, [0.0d0, 1.0d0], [1.0d0])
    solution = solver%get_solution()
    rst = solver%get_order() >= 1 .and. solver%get_order() <= 12 .and. &
        abs(solution(size(solution,1),2) - exp(1.0d0)) < 1.0d-5
end function

function test_bdf_mass_matrix() result(rst)
    logical :: rst
    type(bdf) :: solver
    type(ode_container) :: model
    real(real64), allocatable :: solution(:,:)
    model%fcn => mass_exponential_model
    model%mass_matrix => diagonal_mass_matrix
    call solver%solve(model, [0.0d0, 1.0d0], [1.0d0])
    solution = solver%get_solution()
    rst = abs(solution(size(solution,1),2) - exp(-0.5d0)) < 1.0d-4
end function

subroutine exponential_model(x, y, dydx, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args
    dydx = y
end subroutine

subroutine mass_exponential_model(x, y, dydx, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args
    dydx = -y
end subroutine

subroutine diagonal_mass_matrix(x, y, matrix, args)
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: matrix(:,:)
    class(*), intent(inout), optional :: args
    matrix = 0.0d0
    matrix(1,1) = 2.0d0
end subroutine
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

    if (size(sol, 1) <= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_1 - no internal steps"
    end if

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

    if (size(sol, 1) <= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_1 - no internal steps"
    end if

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

    if (size(sol, 1) <= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bdf_1 - no internal steps"
    end if

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