module diffeq_jacobian_tests
    use iso_fortran_env
    use diffeq_models
    use fortran_test_helper
    use diffeq
    implicit none
contains
! ------------------------------------------------------------------------------
function test_fd_jacobian_1() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    type(ode_container) :: obj
    real(real64) :: x, y(2), dydx(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => vanderpol

    ! Define the conditions at which to evaluate the Jacobian
    x = 0.0d0
    y = 0.0d0
    call vanderpol(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = vanderpol_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_1 1-1"
    end if

    ! Try a new set of coordinates
    call random_number(x)
    call random_number(y)
    call vanderpol(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = vanderpol_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 2
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_1 1-2"
    end if

    ! Format
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_fd_jacobian_2() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    type(ode_container) :: obj
    real(real64) :: x, y(2), dydx(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => duffing

    ! Define the conditions at which to evaluate the Jacobian
    x = 0.0d0
    y = 0.0d0
    call duffing(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = duffing_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_2 1-1"
    end if

    ! Try a new set of coordinates
    call random_number(x)
    call random_number(y)
    call duffing(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = duffing_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 2
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_2 1-2"
    end if

    ! Format
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_fd_jacobian_3() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    type(ode_container) :: obj
    real(real64) :: x, y(2), dydx(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => mathieu

    ! Define the conditions at which to evaluate the Jacobian
    x = 0.0d0
    y = 0.0d0
    call mathieu(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = mathieu_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_3 1-1"
    end if

    ! Try a new set of coordinates
    call random_number(x)
    call random_number(y)
    call mathieu(x, y, dydx)

    ! Compute the solution matrix and the Jacobian estimate
    ans = mathieu_jacobian(x, y)
    call obj%compute_jacobian(x, y, dydx, jac)

    ! Test 2
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_3 1-2"
    end if

    ! Format
100 format(A)
end function

! ------------------------------------------------------------------------------
end module