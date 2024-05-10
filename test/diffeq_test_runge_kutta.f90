module diffeq_test_runge_kutta
    use iso_fortran_env
    use diffeq
    use fortran_test_helper
    use diffeq_models
    implicit none

contains
! ------------------------------------------------------------------------------
function test_runge_kutta_45_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(runge_kutta_45) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_runge_kutta_45_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    type(runge_kutta_45) :: integrator
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
        print "(A)", "TEST FAILED: test_runge_kutta_45_2 -1"
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module