module diffeq_frf_tests
    use iso_fortran_env
    use fortran_test_helper
    use diffeq_harmonics
    use diffeq_models
    implicit none

    type, extends(harmonic_ode_container) :: test_sweep_ode
    contains
        procedure, public :: ode => test_sweep_ode_fcn
    end type

contains
! ------------------------------------------------------------------------------
function test_frf_1() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: tol = 0.05d0
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    type(forced_ode_container) :: mdl
    integer(int32) :: i, ind
    real(real64), allocatable :: freq(:), mag(:), magans(:), ratio(:), &
        reference(:)
    complex(real64), allocatable :: sol(:,:), tf(:), s(:)

    ! Initialization
    rst = .true.
    mdl%fcn => example_2nd_order
    mdl%forcing_function => example_2nd_order_forcing
    
    ! Compute the FRF's
    sol = frequency_response(mdl, 5.0d0, [0.0d0, 0.0d0], freq = freq)
    mag = abs(sol(:,1))

    ! Compute the solution
    s = j * (2.0d0 * pi * freq)
    tf = 1.0d0 / (s**2 + 2.0d0 * z * wn * s + wn**2)
    magans = abs(tf)

    ! The solution will carry on beyond the desired max frequency; however, is
    ! effectively invalid beyond such frequency.  As a result we must limit the
    ! test range
    ind = 0
    do i = 1, size(mag)
        if (freq(i) > fmax) then
            ind = i - 1
            exit
        end if
    end do

    ! Test the magnitude
    ratio = mag(:ind) / magans(:ind)
    allocate(reference(size(ratio)), source = 1.0d0)
    if (.not.assert(ratio, reference, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_frf_1 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
subroutine test_sweep_ode_fcn(this, x, y, dydx)
    ! Arguments
    class(test_sweep_ode), intent(in) :: this
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)

    ! Model Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    real(real64), parameter :: F = 1.0d2
    
    ! Process
    dydx(1) = y(2)
    dydx(2) = F * sin(2.0d0 * pi * this%excitation_frequency * x) - &
        2.0d0 * z * wn * y(2) - wn**2 * y(1)
end subroutine

! ------------------------------------------------------------------------------
function test_frf_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 0.05d0
    integer(int32), parameter :: npts = 100
    real(real64), parameter :: fmin = 1.0d1
    real(real64), parameter :: fmax = 1.0d2
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: z = 1.0d-1
    real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    real(real64), parameter :: F = 1.0d2
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    type(test_sweep_ode) :: sys
    integer(int32) :: i
    real(real64) :: df, freq(npts)
    real(real64), allocatable :: mag(:), magans(:), ratio(:), reference(:)
    complex(real64), allocatable :: sol(:,:), tf(:), s(:)

    ! Initialization
    rst = .true.
    df = (fmax - fmin) / (npts - 1.0d0)
    freq = (/ (df * i + fmin, i = 0, npts - 1) /)
    s = j * (2.0d0 * pi * freq)

    ! Compute the frequency response
    sol = frequency_response(sys, freq, [0.0d0, 0.0d0])
    mag = abs(sol(:,1))

    ! Compute the solution
    tf = F / (s**2 + 2.0d0 * z * wn * s + wn**2)
    magans = abs(tf)

    ! Test the magnitude
    ratio = mag / magans
    allocate(reference(size(ratio)), source = 1.0d0)
    if (.not.assert(ratio, reference, tol)) then
        rst = .false.
        print 100, "TEST FAILED: test_frf_2 1-1"
    end if

    ! Formatting
100 format(A)
end function

! ------------------------------------------------------------------------------
end module