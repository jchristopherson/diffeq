module diffeq_dprk45_tests
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fortran_test_helper

    ! Dormand-Prince 4th/5th Order Model Coefficients
    real(real64), parameter :: a21 = 1.0d0 / 5.0d0
    real(real64), parameter :: a31 = 3.0d0 / 40.0d0
    real(real64), parameter :: a32 = 9.0d0 / 40.0d0
    real(real64), parameter :: a41 = 44.0d0 / 45.0d0
    real(real64), parameter :: a42 = -56.0d0 / 15.0d0
    real(real64), parameter :: a43 = 32.0d0 / 9.0d0
    real(real64), parameter :: a51 = 1.9372d4 / 6.561d3
    real(real64), parameter :: a52 = -2.536d4 / 2.187d3
    real(real64), parameter :: a53 = 6.4448d4 / 6.561d3
    real(real64), parameter :: a54 = -2.12d2 / 7.29d2
    real(real64), parameter :: a61 = 9.017d3 / 3.168d3
    real(real64), parameter :: a62 = -3.55d2 / 33.0d0
    real(real64), parameter :: a63 = 4.6732d4 / 5.247d3
    real(real64), parameter :: a64 = 49.0d0 / 1.76d2
    real(real64), parameter :: a65 = -5.103d3 / 1.8656d4
    real(real64), parameter :: a71 = 35.0d0 / 3.84d2
    real(real64), parameter :: a72 = 0.0d0
    real(real64), parameter :: a73 = 5.0d2 / 1.113d3
    real(real64), parameter :: a74 = 1.25d2 / 1.92d2
    real(real64), parameter :: a75 = -2.187d3 / 6.784d3
    real(real64), parameter :: a76 = 11.0d0 / 84.0d0
    
    real(real64), parameter :: e1 = -71.0d0 / 5.76d4
    real(real64), parameter :: e2 = 0.0d0
    real(real64), parameter :: e3 = 71.0d0 / 1.6695d4
    real(real64), parameter :: e4 = -71.0d0 / 1.92d3
    real(real64), parameter :: e5 = 1.7253d4 / 3.392d5
    real(real64), parameter :: e6 = -22.0d0 / 5.25d2
    real(real64), parameter :: e7 = 1.0d0 / 4.0d1

    real(real64), parameter :: c2 = 1.0d0 / 5.0d0
    real(real64), parameter :: c3 = 3.0d0 / 1.0d1
    real(real64), parameter :: c4 = 4.0d0 / 5.0d0
    real(real64), parameter :: c5 = 8.0d0 / 9.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: c7 = 1.0d0

contains
! ------------------------------------------------------------------------------
function test_dprk45_attempt_step() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: h = 1.0d-1

    ! Local Variables
    type(dprk45_integrator) :: integrator
    type(ode_container) :: mdl
    real(real64) :: x, y(1), yn(1), en(1), yt(1), dydx(1), k2(1), k3(1), &
        k4(1), k5(1), k6(1), yans(1), eans(1), dydxn(1)

    ! Define the model
    mdl%fcn => test_1dof_1

    ! Initialization
    rst = .true.
    call random_number(x)
    call random_number(y)
    call mdl%fcn(x, y, dydx)

    ! Perform the actual steps manually
    yt = y + h * a21 * dydx
    call mdl%fcn(x + c2 * h, yt, k2)

    yt = y + h * (a31 * dydx + a32 * k2)
    call mdl%fcn(x + c3 * h, yt, k3)

    yt = y + h * (a41 * dydx + a42 * k2 + a43 * k3)
    call mdl%fcn(x + c4 * h, yt, k4)

    yt = y + h * (a51 * dydx + a52 * k2 + a53 * k3 + a54 * k4)
    call mdl%fcn(x + c5 * h, yt, k5)

    yt = y + h * (a61 * dydx + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    call mdl%fcn(x + h, yt, k6)

    yans = y + h * (a71 * dydx + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    call mdl%fcn(x + h, yans, dydxn)
    
    eans = h * (e1 * dydx + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * dydxn)

    ! Attempt a single step
    call integrator%attempt_step(mdl, h, x, y, yn, en)

    ! Test 1
    if (.not.assert(yn, yans, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_dprk45_attempt_step 1-1"
    end if

    ! Test 2
    if (.not.assert(en, eans, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_dprk45_attempt_step 1-2"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_dprk45_step() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: xmax = 1.0d0
    real(real64), parameter :: xinit = 0.0d0
    real(real64), parameter :: yinit = 2.0d0
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(dprk45_integrator) :: integrator
    type(ode_container) :: mdl
    real(real64) :: yn(1), yans

    ! Define the model
    mdl%fcn => test_1dof_1

    ! Initialization
    rst = .true.

    ! Take one step using the dprk45 integrator
    call integrator%set_step_size(xmax)
    call integrator%step(mdl, xinit, xmax, [yinit], yn)

    ! Compute the solution
    yans = test_1dof_solution_1(xinit + integrator%get_step_size())
    
    ! Test
    if (.not.assert([yans], [yn], tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_dprk45_step 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_dprk45_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(dprk45_integrator) :: integrator
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
        print 100, "TEST FAILED: test_dprk45_1 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_dprk45_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(dprk45_integrator) :: integrator
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
        print 100, "TEST FAILED: test_dprk45_2 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_dprk45_3() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(dprk45_integrator) :: integrator
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
        print 100, "TEST FAILED: test_dprk45_3 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
end module