module diffeq_equilibrium_tests
    use iso_fortran_env
    use diffeq
    use diffeq_tools
    use diffeq_models
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_equilibrium_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6
    real(real64), parameter :: ans1 = -2.0d0
    real(real64), parameter :: ans2 = 3.0d0

    ! Local Variables
    type(ode_container) :: mdl
    real(real64) :: sol(1,2)

    ! Test
    rst = .true.
    mdl%fcn => first_order_1
    sol = find_equilibrium_points(mdl, [0.0d0])

    ! Results
    if (.not.assert(sol(1,1), ans1, tol) .and. &
        .not.assert(sol(1,1), ans2, tol)) &
    then
        rst = .false.
        print 100, "TEST FAILED: test_equilibrium_1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_equilibrium_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6
    real(real64), parameter :: ans1 = -2.0d0
    real(real64), parameter :: ans2 = -1.0d0
    real(real64), parameter :: ans3 = 2.0d0

    ! Local Variables
    type(ode_container) :: mdl
    real(real64) :: sol(1,2)

    ! Test
    rst = .true.
    mdl%fcn => first_order_2
    sol = find_equilibrium_points(mdl, [-5.0d0])

    ! Results
    if (.not.assert(sol(1,1), ans1, tol) .and. &
        .not.assert(sol(1,1), ans2, tol) .and. &
        .not.assert(sol(1,1), ans3, tol)) &
    then
        rst = .false.
        print 100, "TEST FAILED: test_equilibrium_2"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
end module