module diffeq_bsrk32_tests
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_bsrk32_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(bsrk32_integrator) :: integrator
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
        print 100, "TEST FAILED: test_bsrk32_1 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
function test_bsrk32_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(bsrk32_integrator) :: integrator
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
        print 100, "TEST FAILED: test_bsrk32_2 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
end module