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
    real(real64) :: x, y(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => vanderpol

    ! Define the conditions at which to evaluate the Jacobian
    x = 0.0d0
    y = 0.0d0

    ! Compute the solution matrix and the Jacobian estimate
    ans = vanderpol_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_1 1-1"
    end if

    ! Try a new set of coordinates
    x = 1.5d0
    y = 1.5d0

    ! Compute the solution matrix and the Jacobian estimate
    ans = vanderpol_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

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
    real(real64) :: x, y(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => duffing

    ! Define the conditions at which to evaluate the Jacobian
    x = 0.0d0
    y = 0.0d0

    ! Compute the solution matrix and the Jacobian estimate
    ans = duffing_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_2 1-1"
    end if

    ! Try a new set of coordinates
    x = 1.5d0
    y = 1.5d0

    ! Compute the solution matrix and the Jacobian estimate
    ans = duffing_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

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
    real(real64) :: x, y(2), jac(2, 2), ans(2, 2)

    ! Tolerances
    real(real64), parameter :: tol = 1.0d-6

    ! Initialization
    rst = .true.
    obj%fcn => mathieu

    ! Define the conditions at which to evaluate the Jacobian
    x = 1.5d0
    y = 1.5d0

    ! Compute the solution matrix and the Jacobian estimate
    ans = mathieu_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

    ! Test 1
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_3 1-1"
    end if

    ! Try a new set of coordinates
    call random_number(x)
    call random_number(y)

    ! Compute the solution matrix and the Jacobian estimate
    ans = mathieu_jacobian(x, y)
    call obj%compute_jacobian(x, y, jac)

    ! Test 2
    if (.not.assert(ans, jac, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_3 1-2"
    end if

    ! Format
100 format(A)
end function

! ------------------------------------------------------------------------------
subroutine dummy_jacobian_routine(x, y, jac)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: jac(:,:)

    ! Process
    jac = reshape([1.0d0, 2.0d0, 3.0d0, 4.0d0], [2, 2])
end subroutine

! This test is intended to ensure the user-defined Jacobian routine works as 
! intended, and does not require the definition of an ODE function.
function test_fd_jacobian_4() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: ans(2, 2), jac(2, 2), x, y(2)
    type(ode_container) :: obj

    ! Initialization
    rst = .true.
    x = 0.0d0
    y = 0.0d0
    obj%jacobian => dummy_jacobian_routine

    ! Compute the answer
    call dummy_jacobian_routine(x, y, ans)

    ! Test
    call obj%compute_jacobian(x, y, jac)
    if (.not.assert(ans, jac)) then
        rst = .false.
        print 100, "TEST FAILED: test_fd_jacobian_4 1-1"
    end if

    ! Format
100 format(A)
end function

! ------------------------------------------------------------------------------
end module