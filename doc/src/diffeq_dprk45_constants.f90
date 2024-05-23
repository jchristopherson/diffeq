module diffeq_dprk45_constants
    use iso_fortran_env
    implicit none

    ! Dormand-Prince 4th/5th Order Model Coefficients
    real(real64), parameter :: a21 = 1.0d0 / 5.0d0
    real(real64), parameter :: a31 = 3.0d0 / 40.0d0
    real(real64), parameter :: a32 = 9.0d0 / 40.0d0
    real(real64), parameter :: a41 = 44.0d0 / 45.0d0
    real(real64), parameter :: a42 = -56.0d0 / 15.0d0
    real(real64), parameter :: a43 = 32.0d0 / 9.0d0
    real(real64), parameter :: a51 = 1.9372d4 / 6.561d3
    real(real64), parameter :: a52 = -2.536d4 / 2.187d3
    real(real64), parameter :: a53 = 6.4448d4 / 6.561d3
    real(real64), parameter :: a54 = -2.12d2 / 7.29d2
    real(real64), parameter :: a61 = 9.017d3 / 3.168d3
    real(real64), parameter :: a62 = -3.55d2 / 33.0d0
    real(real64), parameter :: a63 = 4.6732d4 / 5.247d3
    real(real64), parameter :: a64 = 49.0d0 / 1.76d2
    real(real64), parameter :: a65 = -5.103d3 / 1.8656d4
    real(real64), parameter :: a71 = 35.0d0 / 3.84d2
    real(real64), parameter :: a72 = 0.0d0
    real(real64), parameter :: a73 = 5.0d2 / 1.113d3
    real(real64), parameter :: a74 = 1.25d2 / 1.92d2
    real(real64), parameter :: a75 = -2.187d3 / 6.784d3
    real(real64), parameter :: a76 = 11.0d0 / 84.0d0
    
    real(real64), parameter :: e1 = -71.0d0 / 5.76d4
    real(real64), parameter :: e2 = 0.0d0
    real(real64), parameter :: e3 = 71.0d0 / 1.6695d4
    real(real64), parameter :: e4 = -71.0d0 / 1.92d3
    real(real64), parameter :: e5 = 1.7253d4 / 3.392d5
    real(real64), parameter :: e6 = -22.0d0 / 5.25d2
    real(real64), parameter :: e7 = 1.0d0 / 4.0d1

    real(real64), parameter :: c2 = 1.0d0 / 5.0d0
    real(real64), parameter :: c3 = 3.0d0 / 1.0d1
    real(real64), parameter :: c4 = 4.0d0 / 5.0d0
    real(real64), parameter :: c5 = 8.0d0 / 9.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: c7 = 1.0d0

    ! Interpolation Parameters
    real(real64), parameter :: d1 = -1.2715105075d10 / 1.1282082432d10
    real(real64), parameter :: d3 = 8.74874797d10 / 3.2700410799d10
    real(real64), parameter :: d4 = -1.0690763975d10 / 1.880347072d9
    real(real64), parameter :: d5 = 7.01980252875d11 / 1.99316789632d11
    real(real64), parameter :: d6 = -1.453857185d9 / 8.22651844d8
    real(real64), parameter :: d7 = 6.9997945d7 / 2.9380423d7
    
end module