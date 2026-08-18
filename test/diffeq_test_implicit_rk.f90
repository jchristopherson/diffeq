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