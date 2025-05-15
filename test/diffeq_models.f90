module diffeq_models
    use iso_fortran_env
    implicit none
contains
! Van Der Pol Equation
pure subroutine vanderpol(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Model Constants
    real(real64), parameter :: mu = 5.0d0

    ! Equations
    dydx(1) = y(2)
    dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine

function vanderpol_jacobian(x, y) result(rst)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64) :: rst(2, 2)

    ! Model Constants
    real(real64), parameter :: mu = 5.0d0

    ! Jacobian
    rst(1,1) = 0.0d0
    rst(2,1) = -2.0d0 * mu * y(1) * y(2) - 1.0d0
    rst(1,2) = 1.0d0
    rst(2,2) = mu * (1.0d0 - y(1)**2)
end function

subroutine vanderpol_args(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Get the model constant
    real(real64) :: mu
    select type (args)
    type is (real(real64))
        mu = args
    end select

    ! Equations
    dydx(1) = y(2)
    dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine

! ------------------------------------------------------------------------------
! Duffing Equation
pure subroutine duffing(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Model Constants
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 5.0d0
    real(real64), parameter :: delta = 2.0d-2
    real(real64), parameter :: gamma = 8.0d0
    real(real64), parameter :: w = 0.5d0

    ! Equations
    dydx(1) = y(2)
    dydx(2) = gamma * cos(w * x) - delta * y(2) - alpha * y(1) - beta * y(1)**3
end subroutine

function duffing_jacobian(x, y) result(rst)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64) :: rst(2, 2)

    ! Model Constants
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 5.0d0
    real(real64), parameter :: delta = 2.0d-2
    real(real64), parameter :: gamma = 8.0d0
    real(real64), parameter :: w = 0.5d0

    ! Jacobian
    rst(1,1) = 0.0d0
    rst(2,1) = -3.0d0 * beta * y(1)**2 - alpha
    rst(1,2) = 1.0d0
    rst(2,2) = -delta
end function

! ------------------------------------------------------------------------------
pure subroutine mathieu(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Model Constants
    real(real64), parameter :: a = 1.0d0
    real(real64), parameter :: q = 0.2d0

    ! Equations
    dydx(1) = y(2)
    dydx(2) = (2.0d0 * q * cos(2.0d0 * x) - a) * y(1)
end subroutine

function mathieu_jacobian(x, y) result(rst)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64) :: rst(2, 2)

    ! Model Constants
    real(real64), parameter :: a = 1.0d0
    real(real64), parameter :: q = 0.2d0

    ! Jacobian
    rst(1,1) = 0.0d0
    rst(2,1) = 2.0d0 * q * cos(2.0d0 * x) - a
    rst(1,2) = 1.0d0
    rst(2,2) = 0.0d0
end function

! ------------------------------------------------------------------------------
! Linear Test Problem:
! y" + wn**2 * y = 0
! y(0) = 1
! y'(0) = 1/2
! y(x) = sin(wn * x) / 2 / wn + cos(wn * x)
pure subroutine test_2dof_1(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Model Constants
    real(real64), parameter :: wn = 2.0d1

    ! Equations
    dydx(1) = y(2)
    dydx(2) = -wn**2 * y(1)
end subroutine

pure elemental function test_2dof_solution_1(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Model Constants
    real(real64), parameter :: wn = 2.0d1

    ! Solution
    rst = sin(wn * x) / 2.0d0 / wn + cos(wn * x)
end function

! ------------------------------------------------------------------------------
! 1 DOF Test Problem
! y' + y * sin(x)**2 = 0
! y(0) = 2
! y(x) = 2 * exp(0.25 * sin(2 * x) - 0.5 * x)
pure subroutine test_1dof_1(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Equation
    dydx(1) = -y(1) * sin(x)**2
end subroutine

pure elemental function test_1dof_solution_1(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Solution
    rst = 2.0d0 * exp(0.25d0 * sin(2.0d0 * x) - 0.5d0 * x)
end function

! ------------------------------------------------------------------------------
! 2nd Order Test Problem
! x" + 2 * z * wn * x' + wn**2 * x = f(t)
! pure function example_2nd_order_forcing(t) result(rst)
!     use diffeq_harmonics, only : chirp

!     ! Arguments
!     real(real64), intent(in) :: t
!     real(real64) :: rst

!     ! Process
!     rst = chirp(t, 1.0d2, 5.0d0, 1.0d0, 1.0d2)
! end function

! subroutine example_2nd_order(t, x, dxdt)
!     ! Arguments
!     real(real64), intent(in) :: t
!     real(real64), intent(in), dimension(:) :: x
!     real(real64), intent(out), dimension(:) :: dxdt

!     ! Model Parameters
!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
!     real(real64), parameter :: z = 1.0d-1
!     real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1

!     ! Local Variables
!     real(real64) :: f

!     ! Process
!     f = example_2nd_order_forcing(t)
!     dxdt(1) = x(2)
!     dxdt(2) = f - (2.0d0 * z * wn * x(2) + wn**2 * x(1))
! end subroutine

! ------------------------------------------------------------------------------
! Roots: -2, 3
pure subroutine first_order_1(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    ! Process
    dydx(1) = y(1)**2 - y(1) - 6.0d0
end subroutine

! ------------------------------------------------------------------------------
! Roots: -2, 2, -1
pure subroutine first_order_2(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx
    class(*), intent(inout), optional :: args

    ! Process
    dydx(1) = (y(1)**2 - 4.0d0) * (y(1) + 1.0d0)**2
end subroutine

! ------------------------------------------------------------------------------
subroutine parametric_model(t, x, dxdt, args)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt
    class(*), intent(inout), optional :: args

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: m = 6.1763013134581447d-5
    real(real64), parameter :: b = 4.3487216190266362d-3
    real(real64), parameter :: k = 190.46442027195485d0
    real(real64), parameter :: kappa = 3.5035644222521636d-3
    real(real64), parameter :: d = 8.0d0
    real(real64), parameter :: V = 3.0d1
    real(real64), parameter :: startRatio = 5.0d-1
    real(real64), parameter :: endRatio = 2.0d0
    real(real64), parameter :: span = 1.0d0

    ! Local Variables
    real(real64) :: wn, fn, startfHz, endfHz, c, F, y

    ! Forcing Routine
    wn = sqrt(k / m)
    fn = wn / (2.0d0 * pi)
    startfHz = startRatio * fn
    endfHz = endRatio * fn
    c = (endfHz - startfHz) / span
    y = V * sin(2.0d0 * pi * t * (0.5d0 * c * t + startfHz))
    F = sign(1.0d0, y) * kappa * y**2 / ((1.0d0 - x(1) / d)**2)

    ! Equations of motion:
    ! m x" + b x' + k x = sign(y) kappa y**2 / (1 - x / d)**2
    dxdt(1) = x(2)
    dxdt(2) = (F - b * x(2) - k * x(1)) / m
end subroutine

! ------------------------------------------------------------------------------
end module