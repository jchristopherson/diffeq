submodule (diffeq_harmonics) diffeq_harmonics_routines
    use fftpack
    use spectrum
contains
! ------------------------------------------------------------------------------
module function frf_fft_1(sys, span, iv, fs, solver, win, freq, &
    err) result(rst)
    ! Arguments
    class(forced_ode_container), intent(inout) :: sys
    real(real64), intent(in) :: span
    real(real64), intent(in), dimension(:) :: iv
    real(real64), intent(in), optional :: fs
    class(ode_integrator), intent(inout), optional, target :: solver
    class(window), intent(in), optional, target :: win
    real(real64), intent(out), optional, allocatable, dimension(:) :: freq
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: i, npts, neqn, nfreq, flag
    real(real64) :: dt, samplerate, df
    real(real64), allocatable, dimension(:) :: t, force
    real(real64), allocatable, dimension(:,:) :: sol
    class(ode_integrator), pointer :: integrator
    type(dprk45_integrator), target :: default_integrator
    class(window), pointer :: w
    type(hamming_window), target :: default_window
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(fs)) then
        samplerate = fs
    else
        samplerate = 1.024d3
    end if
    if (present(solver)) then
        integrator => solver
    else
        integrator => default_integrator
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    dt = 1.0d0 / samplerate
    npts = floor(span / dt) + 1
    neqn = size(iv)
    if (present(win)) then
        w => win
    else
        w => default_window
        w%size = npts
    end if
    nfreq = compute_transform_length(w%size)

    ! Input Checking
    if (.not.associated(sys%forcing_function)) go to 20

    ! Memory Allocation
    allocate(rst(nfreq, neqn), stat = flag, source = (0.0d0, 0.0d0))
    if (flag == 0) allocate(t(npts), stat = flag)
    if (flag == 0) allocate(force(npts), stat = flag)
    if (flag == 0 .and. present(freq)) allocate(freq(nfreq), stat = flag)
    if (flag /= 0) go to 10

    ! Construct the time array and forcing function
    do i = 1, npts
        t(i) = dt * (i - 1.0d0)
        force(i) = sys%forcing_function(t(i))
    end do

    ! Solve the ODE's
    sol = integrator%solve(sys, t, iv, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the transfer functions of each result
    do i = 1, neqn
        rst(:,i) = siso_transfer_function(w, force, sol(:,i+1), &
            err = errmgr)
        if (errmgr%has_error_occurred()) return
    end do

    ! See if we need to build a frequency vector
    if (present(freq)) then
        df = frequency_bin_width(samplerate, w%size)
        freq = (/ (df * i, i = 0, nfreq - 1) /)
    end if

    ! End
    return

    ! Memory Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("frf_fft_1", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! No forcing function defined
20  continue
    call errmgr%report_error("frf_fft_1", "No forcing routine defined.", &
        DIFFEQ_NULL_POINTER_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
module function frf_sweep(sys, freq, iv, solver, ncycles, ntransient, &
    pointspercycle, err) result(rst)
    ! Arguments
    class(harmonic_ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: freq, iv
    class(ode_integrator), intent(inout), optional, target :: solver
    integer(int32), intent(in), optional :: ncycles, ntransient, pointspercycle
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable, dimension(:,:) :: rst

    ! Parameters
    real(real64), parameter :: zerotol = sqrt(epsilon(0.0d0))

    ! Local Variables
    integer(int32) :: i, j, nfreq, neqn, nc, nt, ntotal, npts, ppc, flag, &
        nfft, i1, ncpts
    real(real64) :: dt, tare, phase, amp
    real(real64), allocatable, dimension(:) :: ic, t
    real(real64), allocatable, dimension(:,:) :: sol
    complex(real64), allocatable, dimension(:) :: xpts
    class(ode_integrator), pointer :: integrator
    type(dprk45_integrator), target :: default_integrator
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(ncycles)) then
        nc = ncycles
    else
        nc = 20
    end if
    if (present(ntransient)) then
        nt = ntransient
    else
        nt = 200
    end if
    if (present(pointspercycle)) then
        ppc = pointspercycle
    else
        ppc = 1000
    end if
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nfreq = size(freq)
    neqn = size(iv)
    ntotal = nt + nc
    npts = ntotal * ppc
    ncpts = nc * ppc
    i1 = npts - ncpts + 1
    nfft = 2**next_power_of_two(ppc * nc)

    ! Set up the integrator
    if (present(solver)) then
        integrator => solver
    else
        integrator => default_integrator
    end if

    ! Input Checking
    if (nc < 1) go to 20
    if (nt < 1) go to 30
    if (ppc < 2) go to 40
    do i = 1, nfreq
        if (abs(freq(i)) < zerotol) go to 50
    end do

    ! Local Memory Allocation
    allocate(rst(nfreq, neqn), stat = flag)
    if (flag == 0) allocate(ic(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(t(npts), stat = flag)
    if (flag == 0) allocate(xpts(nfft), stat = flag, source = (0.0d0, 0.0d0))
    if (flag /= 0) go to 10

    ! Cycle over each frequency point
    do i = 1, nfreq
        ! Define the time vector
        dt = (1.0d0 / freq(i)) / (ppc - 1.0d0)
        t = (/ (dt * j, j = 0, npts - 1) /)

        ! Update the frequency
        sys%excitation_frequency = freq(i)

        ! Compute the solution
        sol = integrator%solve(sys, t, ic, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Reset the initial conditions to the last solution point
        ic = sol(npts, 2:)

        ! Determine the magnitude and phase for each equation
        do j = 1, neqn
            rst(i,j) = get_magnitude_phase(sol(i1:,j+1), xpts)
        end do
    end do

    ! Offset phase to be zero at the lowest frequency for each DOF
    i1 = minloc(freq, 1)
    do j = 1, neqn
        tare = atan2(aimag(rst(i1,j)), real(rst(i1,j)))
        do i = 1, nfreq
            phase = atan2(aimag(rst(i,j)), real(rst(i,j))) - tare
            amp = abs(rst(i,j))
            rst(i,j) = cmplx(amp * cos(phase), amp * sin(phase))
        end do
    end do

    ! End
    return

    ! Memory Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("frf_sweep", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Number of Cycles Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The number of cycles to analyze must be at " // &
        "least 1; however, a value of ", nc, " was found."
    call errmgr%report_error("frf_sweep", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Number of Transient Cycles Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The number of transient cycles must be at " // &
        "least 1; however, a value of ", nt, " was found."
    call errmgr%report_error("frf_sweep", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Points Per Cycle Error
40  continue
    write(errmsg, 100) "The number of points per cycle must be at " // &
        "least 2; however, a value of ", ppc, " was found."
    call errmgr%report_error("frf_sweep", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Zero-Valued Frequency Error
50  continue
    write(errmsg, 100) "A zero-valued frequency was found at index ", i, "."
    call errmgr%report_error("frf_sweep", trim(errmsg), &
        DIFFEQ_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
function get_magnitude_phase(x, xzeros) result(rst)
    ! Arguments
    real(real64), intent(in), dimension(:) :: x
    complex(real64), intent(inout), dimension(:) :: xzeros
    complex(real64) :: rst

    ! Local Variables
    integer(int32) :: ind, n, m, nx
    real(real64) :: amp, phase

    ! Initialization
    nx = size(x)
    n = size(xzeros)
    m = compute_transform_length(n)
    amp = 0.5d0 * (maxval(x) - minval(x))

    ! Zero pad the data
    xzeros(1:nx) = cmplx(x, 0.0d0, real64)
    xzeros(nx+1:) = cmplx(0.0d0, 0.0d0, real64)

    ! Compute the FFT to estimate the phase
    xzeros = fft(xzeros)
    ind = maxloc(abs(xzeros(1:m)), 1)
    phase = atan2(aimag(xzeros(ind)), real(xzeros(ind)))
    rst = cmplx(amp * cos(phase), amp * sin(phase), real64)
end function

! ------------------------------------------------------------------------------
end submodule