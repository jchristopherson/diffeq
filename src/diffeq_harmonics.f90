!> @brief This module contains routines for assisting in the harmonic analysis
!! of systems of ODE's.
module diffeq_harmonics
    use iso_fortran_env
    use diffeq
    use ferror
    use spectrum
    implicit none
    private
    public :: ode_excite
    public :: forced_ode_container
    public :: harmonic_ode_container
    public :: chirp
    public :: frequency_response

    interface
        function ode_excite(t) result(rst)
            use iso_fortran_env
            real(real64), intent(in) :: t
            real(real64) :: rst
        end function
    end interface

    type, extends(ode_container) :: forced_ode_container
        procedure(ode_excite), pointer, nopass :: forcing_function
    end type

    type, extends(ode_container) :: harmonic_ode_container
        real(real64) :: excitation_frequency
    end type

    !> @brief Computes the frequency response of each equation in a system of
    !! forced differential equations.
    !!
    !! @par Syntax 1
    !! This implementation utilizes a discrete Fourier transform to estimate the
    !! frequency response functions for each equation in a system of ODE's.
    !! @code{.f90}
    !! allocatable complex(real64)(:,:) function frequency_response ( &
    !!  class(forced_ode_container) sys, &
    !!  real(real64) span, &
    !!  real(real64) iv(:), &
    !!  optional real(real64) fs, &
    !!  optional class(ode_integrator) solver, &
    !!  optional class(window) win, &
    !!  optional allocatable real(real64) freq(:), &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in,out] sys The @p forced_ode_container containing the system of
    !!  differential equations and associated forcing function to analyze.
    !! @param[in] span The duration of the time analysis.
    !! @param[in] iv An N-element array containing the initial conditions for
    !!  each of the N equations.
    !! @param[in] fs An optional rate at which to sample the differential
    !!  equation solution.  The default is 1024 Hz, assuming that @p span is
    !!  provided in units of seconds.
    !! @param[in,out] solver An optional differential equation solver.  The 
    !!  default solver is the Dormand-Prince Runge-Kutta integrator
    !!  @ref dprk45_integrator.
    !! @param[in] win An optional parameter allowing the differential equation
    !!  solution to be windowed as part of the FFT process used to construct the
    !!  system transfer functions.  The default is a Hamming window sized to
    !!  contain the entire solution in one window.
    !! @param[out] freq An optional allocatable array that, if supplied, 
    !!  will contain the frequency values at which the solution is provided.
    !!  The units will match the units of @p fs.
    !! @param[in,out] err An optional errors-based object that if provided 
    !!  can be used to retrieve information relating to any errors 
    !!  encountered during execution. If not provided, a default 
    !!  implementation of the errors class is used internally to provide 
    !!  error handling. Possible errors and warning messages that may be 
    !!  encountered are as follows.
    !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
    !!      allocation issue.
    !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined, or
    !!      if no forcing function is defined.
    !!
    !! @return An M-by-N matrix containing the complex-valued frequency response
    !!  functions normalized by the magnitude of the forcing function for each
    !!  of the N equations in the system of ODEs.
    !!
    !! @par Example
    !! The following example illustrates how to compute the frequency response
    !! of a single-degree-of-freedom mechanical vibrating system with a 50 Hz
    !! resonant frequency.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq_harmonics
    !!     use diffeq_models
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     real(real64), parameter :: fmin = 1.0d0
    !!     real(real64), parameter :: fmax = 1.0d2
    !!     real(real64), parameter :: fs = 2.56d2
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!     real(real64), parameter :: z = 1.0d-1
    !!     real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    !!     complex(real64), parameter :: j = (0.0d0, 1.0d0)
    !!
    !!     ! Local Variables
    !!     type(forced_ode_container) :: mdl
    !!     real(real64), allocatable :: freq(:), mag(:), phase(:), amag(:), &
    !!         aphase(:)
    !!     complex(real64), allocatable :: sol(:,:), tf(:), s(:)
    !!
    !!     ! Plot Variables
    !!     type(multiplot) :: plt
    !!     type(plot_2d) :: plt1, plt2
    !!     type(plot_data_2d) :: pd1, pd2
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Define the model
    !!     mdl%fcn => example_2nd_order
    !!     mdl%forcing_function => example_2nd_order_forcing
    !!
    !!     ! Compute the FRF's
    !!     sol = frequency_response(mdl, 5.0d0, [0.0d0, 0.0d0], freq = freq)
    !!
    !!     ! Compute the magnitude and phase
    !!     mag = 2.0d1 * log10(abs(sol(:,1) / sol(1,1)))
    !!     phase = (1.8d2 / pi) * atan2(aimag(sol(:,1)), real(sol(:,1)))
    !!
    !!     ! Compute the analytical solution
    !!     s = j * (2.0d0 * pi * freq)
    !!     tf = 1.0d0 / (s**2 + 2.0d0 * z * wn * s + wn**2)
    !!     amag = 2.0d1 * log10(abs(tf / tf(1)))
    !!     aphase = (1.8d2 / pi) * atan2(aimag(tf), real(tf))
    !!
    !!     ! Plot the results
    !!     call plt%initialize(2, 1)
    !!     call plt1%initialize()
    !!     xAxis => plt1%get_x_axis()
    !!     yAxis => plt1%get_y_axis()
    !!     lgnd => plt1%get_legend()
    !!     call xAxis%set_title("f [Hz]")
    !!     call yAxis%set_title("|X / F| [dB]")
    !!     call xAxis%set_autoscale(.false.)
    !!     call xAxis%set_limits(0.0d0, 1.0d2)
    !!     call lgnd%set_is_visible(.true.)
    !!
    !!     call pd1%define_data(freq, mag)
    !!     call pd1%set_line_width(2.0)
    !!     call pd1%set_name("Numerical")
    !!     call plt1%push(pd1)
    !!
    !!     call pd2%define_data(freq, amag)
    !!     call pd2%set_line_width(3.0)
    !!     call pd2%set_line_style(LINE_DASHED)
    !!     call pd2%set_name("Analytical")
    !!     call plt1%push(pd2)
    !!
    !!     call plt2%initialize()
    !!     xAxis => plt2%get_x_axis()
    !!     yAxis => plt2%get_y_axis()
    !!     call xAxis%set_title("f [Hz]")
    !!     call yAxis%set_title("{/Symbol f} [rad]")
    !!     call xAxis%set_autoscale(.false.)
    !!     call xAxis%set_limits(0.0d0, 1.0d2)
    !!
    !!     call pd1%define_data(freq, phase)
    !!     call plt2%push(pd1)
    !!
    !!     call pd2%define_data(freq, aphase)
    !!     call plt2%push(pd2)
    !!
    !!     call plt%set(1, 1, plt1)
    !!     call plt%set(2, 1, plt2)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The ODE routine and forcing function routine were stored in a seperate
    !! module; however, here is the code for these routines.
    !! @code{.f90}
    !! pure function example_2nd_order_forcing(t) result(rst)
    !!     use diffeq_harmonics, only : chirp
    !!
    !!     ! Arguments
    !!     real(real64), intent(in) :: t
    !!     real(real64) :: rst
    !!
    !!     ! Process
    !!     rst = chirp(t, 1.0d2, 5.0d0, 1.0d0, 1.0d2)
    !! end function
    !!
    !! subroutine example_2nd_order(t, x, dxdt)
    !!     ! Arguments
    !!     real(real64), intent(in) :: t
    !!     real(real64), intent(in), dimension(:) :: x
    !!     real(real64), intent(out), dimension(:) :: dxdt
    !!
    !!     ! Model Parameters
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!     real(real64), parameter :: z = 1.0d-1
    !!     real(real64), parameter :: wn = 2.0d0 * pi * 5.0d1
    !!
    !!     ! Local Variables
    !!     real(real64) :: f
    !!
    !!     ! Process
    !!     f = example_2nd_order_forcing(t)
    !!     dxdt(1) = x(2)
    !!     dxdt(2) = f - (2.0d0 * z * wn * x(2) + wn**2 * x(1))
    !! end subroutine
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html frf_example_1.png
    interface frequency_response
        module procedure :: frf_fft_1
        module procedure :: frf_sweep
    end interface

    ! diffeq_harmonics_routines.f90
    interface
        module function frf_fft_1(sys, span, iv, fs, solver, win, freq, err) &
            result(rst)
            ! Arguments
            class(forced_ode_container), intent(inout) :: sys
            real(real64), intent(in) :: span
            real(real64), intent(in), dimension(:) :: iv
            real(real64), intent(in), optional :: fs
            class(ode_integrator), intent(inout), optional, target :: solver
            class(window), intent(in), optional, target :: win
            real(real64), intent(out), optional, allocatable, dimension(:) :: &
                freq
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable, dimension(:,:) :: rst
        end function

        module function frf_sweep(sys, freq, iv, solver, ncycles, ntransient, &
            pointspercycle, err) result(rst)
            class(harmonic_ode_container) :: sys
            real(real64), intent(in), dimension(:) :: freq, iv
            class(ode_integrator), intent(inout), optional, target :: solver
            integer(int32), intent(in), optional :: ncycles, ntransient, &
                pointspercycle
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable, dimension(:,:) :: rst
        end function
    end interface

contains
    !> @brief Computes a value of a linear chirp signal.
    !!
    !! @param[in] t The value of the independent variable at which to evaluate
    !!  the chirp.
    !! @param[in] amp The amplitude of the signal.
    !! @param[in] span The duration of time it takes to sweep from @p f1Hz to
    !!  @p f2Hz.
    !! @param[in] f1Hz The lower excitation frequency, in Hz.
    !! @param[in] f2Hz The upper excitation frequency, in Hz.
    !!
    !! @return The resulting value of the function at @p t.
    pure module elemental function chirp(t, amp, span, f1Hz, f2Hz) result(rst)
        ! Arguments
        real(real64), intent(in) :: t, amp, span, f1Hz, f2Hz
        real(real64) :: rst

        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        real(real64) :: c

        ! Process
        c = (f2Hz - f1Hz) / span
        rst = amp * sin(2.0d0 * pi * t * (0.5d0 * c * t + f1Hz))
    end function
end module