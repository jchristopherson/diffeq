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

end module