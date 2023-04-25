module diffeq_rkfixed_tests
    use iso_fortran_env
    use diffeq
    use diffeq_models
    implicit none
contains
! ------------------------------------------------------------------------------
function test_rk4_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-4

    ! Local Variables
    type(rk4_fixed_integrator) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts), ans(npts)

    ! Define the model
    mdl%fcn => test_1dof_1

    ! Define the values of x where the solution should be computed
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Compute the actual solution
    ans = test_1dof_solution_1(x)
end function

! ------------------------------------------------------------------------------
end module