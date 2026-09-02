module diffeq_test_pece
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
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_1dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_1 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_dense() result(rst)
    ! Tests the dense output capability.

    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: h = 1.0d-3
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts)
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)
    x = (/ (i * h, i = 0, npts - 1) /)

    ! Compute the solution
    call integrator%solve(mdl, x, [1.0d0, 0.5d0])
    sol = integrator%get_solution()

    ! The solution must be returned at precisely the requested points
    if (size(sol, 1) /= npts) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_dense -1"
        return
    end if
    if (.not.assert(x, sol(:,1), 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_dense -2"
    end if

    ! Compute the actual solution
    ans = test_2dof_solution_1(sol(:,1))

    ! Test
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_dense -3"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_all_steps() result(rst)
    ! Tests that supplying only the end points of the integration range causes
    ! every successful step to be returned.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-5

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    integer(int32) :: i, n
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()
    n = size(sol, 1)

    ! More than just the end points must be reported
    if (n <= 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_all_steps -1"
        return
    end if

    ! The steps must advance monotonically and span the requested range
    do i = 2, n
        if (sol(i,1) <= sol(i-1,1)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_adams_all_steps -2"
            return
        end if
    end do
    if (abs(sol(1,1)) > 1.0d-12 .or. abs(sol(n,1) - 1.0d0) > 1.0d-12) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_all_steps -3"
    end if

    ! Every reported step must be accurate
    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_all_steps -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_mass_matrix() result(rst)
    ! Tests a system carrying a nonsingular mass matrix.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)
    integer(int32) :: n

    ! Initialization
    rst = .true.
    mdl%fcn => pece_mass_ode
    mdl%mass_matrix => pece_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)
    call integrator%set_absolute_tolerance(1.0d-10)
    call integrator%set_relative_tolerance(1.0d-10)

    ! Compute the solution
    call integrator%solve(mdl, [0.0d0, 1.0d0], [1.0d0, 1.0d0])
    sol = integrator%get_solution()
    n = size(sol, 1)

    ! M * y' = g(y) with M = diag(2,3) reduces to y' = [-y1, -2*y2]
    if (abs(sol(n,2) - exp(-1.0d0)) > tol .or. &
        abs(sol(n,3) - exp(-2.0d0)) > tol) &
    then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_mass_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
subroutine pece_mass_matrix(x, y, m, args)
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
subroutine pece_mass_ode(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    ! M * y' = g(y), with M = diag(2,3) and g(y) = [-2*y1, -6*y2]
    dydx(1) = -2.0d0 * y(1)
    dydx(2) = -6.0d0 * y(2)
end subroutine

! ------------------------------------------------------------------------------
function test_adams_with_args() result(rst)
    ! Tests that user-defined arguments reach the model, and in doing so that
    ! the same integrator object may be reused across models.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl, refmdl
    real(real64) :: mu
    real(real64), allocatable :: sol(:,:), ref(:,:)
    integer(int32) :: i
    real(real64) :: x(50)

    ! Initialization
    rst = .true.

    ! The Van der Pol equation is stiff for large mu, so a modest value is
    ! used here; a non-stiff method is the wrong tool otherwise.
    mu = 1.0d0
    x = (/ (i * 1.0d-1, i = 0, 49) /)
    mdl%fcn => vanderpol_args
    refmdl%fcn => vanderpol_mu1
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)

    ! Solve with the parameter supplied through args
    call integrator%solve(mdl, x, [2.0d0, 0.0d0], mu)
    sol = integrator%get_solution()

    ! Solve the equivalent model carrying the parameter internally
    call integrator%clear_buffer()
    call integrator%solve(refmdl, x, [2.0d0, 0.0d0])
    ref = integrator%get_solution()

    ! Test
    if (.not.assert(sol, ref, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_with_args -1"
    end if
end function

! ------------------------------------------------------------------------------
subroutine vanderpol_mu1(x, y, dydx, args)
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    dydx(1) = y(2)
    dydx(2) = 1.0d0 * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine

! ------------------------------------------------------------------------------
function test_adams_order_range() result(rst)
    ! Tests that the order remains within the range the method supports and
    ! that the order limit may be constrained by the caller.

    ! Arguments
    logical :: rst

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1

    ! The method must report its documented order limit
    if (integrator%get_maximum_supported_order() /= 12) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_order_range -1"
    end if

    ! That limit must also be the default, absent any request from the caller
    if (integrator%get_maximum_order() /= 12) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_order_range -2"
    end if

    ! A request beyond that limit must be clamped
    call integrator%set_maximum_order(15)
    if (integrator%get_maximum_order() /= 12) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_order_range -3"
    end if

    ! Restricting the order must still produce a usable solution
    call integrator%set_maximum_order(2)
    call integrator%set_absolute_tolerance(1.0d-8)
    call integrator%set_relative_tolerance(1.0d-8)
    call integrator%solve(mdl, [0.0d0, 1.0d0], [2.0d0])
    sol = integrator%get_solution()

    if (integrator%get_order() > 2) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_order_range -4"
    end if

    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), 1.0d-4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_order_range -5"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_high_order() result(rst)
    ! A smooth problem held to a tight tolerance must actually drive the order
    ! above what the method previously permitted, and must remain accurate
    ! while doing so.

    ! Arguments
    logical :: rst

    ! Local Variables
    type(adams) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:), ans(:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_1dof_1
    call integrator%set_absolute_tolerance(1.0d-12)
    call integrator%set_relative_tolerance(1.0d-12)

    ! Perform the integration
    call integrator%solve(mdl, [0.0d0, 1.0d1], [2.0d0])
    sol = integrator%get_solution()

    ! The order must have climbed beyond the five the method once allowed
    if (integrator%get_order() <= 5) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_high_order -1"
    end if

    ! ... and the answer must still be right
    ans = test_1dof_solution_1(sol(:,1))
    if (.not.assert(ans, sol(:,2), 1.0d-8)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_high_order -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_adams_matches_bdf() result(rst)
    ! Two integrators of altogether different construction, held to the same
    ! tolerances, must agree upon the solution of a non-stiff problem.

    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: npts = 200
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    type(adams) :: pc
    type(bdf) :: implicit
    type(ode_container) :: mdl
    integer(int32) :: i
    real(real64) :: x(npts)
    real(real64), allocatable :: a(:,:), b(:,:)

    ! Initialization
    rst = .true.
    mdl%fcn => test_2dof_1
    x = (/ (i * 5.0d-3, i = 0, npts - 1) /)

    call pc%set_absolute_tolerance(1.0d-10)
    call pc%set_relative_tolerance(1.0d-10)
    call implicit%set_absolute_tolerance(1.0d-10)
    call implicit%set_relative_tolerance(1.0d-10)

    ! Compute both solutions
    call pc%solve(mdl, x, [1.0d0, 0.5d0])
    a = pc%get_solution()
    call implicit%solve(mdl, x, [1.0d0, 0.5d0])
    b = implicit%get_solution()

    ! Test
    if (.not.assert(a, b, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_adams_matches_bdf -1"
    end if
end function

! ------------------------------------------------------------------------------
end module
