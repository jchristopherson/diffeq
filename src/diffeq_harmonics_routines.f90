submodule (diffeq_harmonics) diffeq_harmonics_routines
    use fftpack
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

    ! Local Variables
    integer(int32) :: i, j, nfreq, neqn, nc, nt, ntotal, npts, ppc, flag, lwork
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: ic, t, fftwork
    real(real64), allocatable, dimension(:,:) :: sol
    class(ode_integrator), pointer :: integrator
    type(dprk45_integrator), target :: default_integrator
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(solver)) then
        integrator => solver
    else
        integrator => default_integrator
    end if
    if (present(ncycles)) then
        nc = ncycles
    else
        nc = 5
    end if
    if (present(ntransient)) then
        nt = ntransient
    else
        nt = 30
    end if
    if (present(pointspercycle)) then
        ppc = pointspercycle
    else
        ppc = 30
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
    lwork = 2 * nc * ppc + 15

    ! Input Checking

    ! Local Memory Allocation
    allocate(rst(nfreq, neqn), stat = flag)
    if (flag == 0) allocate(ic(neqn), stat = flag, source = iv)
    if (flag == 0) allocate(t(npts), stat = flag)
    if (flag == 0) allocate(fftwork(lwork), stat = flag)
    if (flag /= 0) go to 10

    ! Set up the Fourier transform
    call dffti(nc * ppc, fftwork)

    ! Cycle over each frequency point
    do i = 1, nfreq
        ! Define the time vector
        dt = ntotal / freq(i)
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
            rst(i,j) = get_magnitude_phase(sol(nt*ppc+1:,j+1), fftwork)
        end do
    end do

    ! End
    return

    ! Memory Error
10  continue
    return
end function

! ------------------------------------------------------------------------------
function get_magnitude_phase(x, work) result(rst)
    ! Arguments
    real(real64), intent(inout), dimension(:) :: x, work
    complex(real64) :: rst

    ! Local Variables
    integer(int32) :: i, n, m
    real(real64) :: mag, check, scale, w, sumw
    complex(real64) :: val
    type(hamming_window) :: win

    ! Initialization
    n = size(x)
    win%size = n
    m = compute_transform_length(n)
    if (mod(n, 2) == 0) then
        scale = 2.0d0 / n
    else
        scale = 2.0d0 / (n - 1.0d0)
    end if

    ! Window the data
    sumw = 0.0d0
    do i = 1, n
        w = win%evaluate(i - 1)
        x(i) = w * x(i)
        sumw = sumw + w
    end do
    scale = scale * n / sumw

    ! Transform the data and find the largest magnitude component
    call dfftf(n, x, work)
    check = 0.0d0
    do i = 2, m
        val = scale * cmplx(x(2*i-1), x(2*i), real64)
        mag = abs(val)
        if (mag > check) then
            rst = val
            check = mag
        end if
    end do
end function

! ------------------------------------------------------------------------------
end submodule