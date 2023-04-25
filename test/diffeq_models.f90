module diffeq_models
    use iso_fortran_env
    implicit none
contains
! Van Der Pol Equation
subroutine vanderpol(x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

    ! Model Constants
    real(real64), parameter :: mu = 2.0d0

    ! Equations
    dydx(1) = y(2)
    dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine

function vanderpol_jacobian(x, y) result(rst)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64) :: rst(2, 2)

    ! Model Constants
    real(real64), parameter :: mu = 2.0d0

    ! Jacobian
    rst(1,1) = 0.0d0
    rst(2,1) = -2.0d0 * mu * y(1) * y(2) - 1.0d0
    rst(1,2) = 1.0d0
    rst(2,2) = mu * (1.0d0 - y(1)**2)
end function

! ------------------------------------------------------------------------------
! Duffing Equation
subroutine duffing(x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

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
subroutine mathieu(x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

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
subroutine test_2dof_1(x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

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
subroutine test_1dof_1(x, y, dydx)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

    ! Equation
    dydx(1) = -y(1) * sin(2.0d0 * x)
end subroutine

pure elemental function test_1dof_solution_1(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    real(real64) :: rst

    ! Solution
    rst = 2.0d0 * exp(0.25d0 * sin(2.0d0 * x) - 0.5d0 * x)
end function

end module