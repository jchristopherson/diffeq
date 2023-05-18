submodule (diffeq) diffeq_rosenbrock
    use linalg

    ! Model Constants
    real(real64), parameter :: c2 = 0.386d0
    real(real64), parameter :: c3 = 0.21d0
    real(real64), parameter :: c4 = 0.63d0
    real(real64), parameter :: bet2p = 0.0317d0
    real(real64), parameter :: bet3p = 0.0635d0
    real(real64), parameter :: bet4p = 0.3438d0
    real(real64), parameter :: d1 = 0.25d0
    real(real64), parameter :: d2 = -0.1043d0
    real(real64), parameter :: d3 = 0.1035d0
    real(real64), parameter :: d4 = -0.3620000000000023d-1
    real(real64), parameter :: a21 = 0.1544d1
    real(real64), parameter :: a31 = 0.9466785280815826d0
    real(real64), parameter :: a32 = 0.2557011698983284d0
    real(real64), parameter :: a41 = 0.3314825187068521d1
    real(real64), parameter :: a42 = 0.2896124015972201d1
    real(real64), parameter :: a43 = 0.9986419139977817d0
    real(real64), parameter :: a51 = 0.1221224509226641d1
    real(real64), parameter :: a52 = 0.6019134481288629d1
    real(real64), parameter :: a53 = 0.1253708332932087d2
    real(real64), parameter :: a54 = -0.6878860361058950d0
    real(real64), parameter :: c21 = -0.56688d1
    real(real64), parameter :: c31 = -0.2430093356833875d1
    real(real64), parameter :: c32 = -0.2063599157091915d0
    real(real64), parameter :: c41 = -0.1073529058151375d0
    real(real64), parameter :: c42 = -0.9594562251023355d1
    real(real64), parameter :: c43 = -0.2047028614809616d2
    real(real64), parameter :: c51 = 0.7496443313967647d1
    real(real64), parameter :: c52 = -0.1024680431464352d2
    real(real64), parameter :: c53 = -0.3399990352819905d2
    real(real64), parameter :: c54 = 0.1170890893206160d2
    real(real64), parameter :: c61 = 0.8083246795921522d1
    real(real64), parameter :: c62 = -0.7981132988064893d1
    real(real64), parameter :: c63 = -0.3152159432874371d2
    real(real64), parameter :: c64 = 0.1631930543123136d2
    real(real64), parameter :: c65 = -0.6058818238834054d1
    real(real64), parameter :: gamma = 0.25d0
    real(real64), parameter :: d21 = 0.1012623508344586d2
    real(real64), parameter :: d22 = -0.7487995877610167d1
    real(real64), parameter :: d23 = -0.3480091861555747d2
    real(real64), parameter :: d24 = -0.7992771707568823d1
    real(real64), parameter :: d25 = 0.1025137723295662d1
    real(real64), parameter :: d31 = -0.6762803392801253d0
    real(real64), parameter :: d32 = 0.6087714651680015d1
    real(real64), parameter :: d33 = 0.1643084320892478d2
    real(real64), parameter :: d34 = 0.2476722511418386d2
    real(real64), parameter :: d35 = -0.6594389125716872d1

contains
! ------------------------------------------------------------------------------
module subroutine rbk_alloc_workspace(this, neqn, err)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: m, n, flag
    character(len = :), allocatable :: errmsg

    ! Process
    m = neqn
    n = 6
    if (allocated(this%m_work)) then
        if (size(this%m_work, 1) /= m .or. size(this%m_work, 2) /= n) then
            deallocate(this%m_work)
            allocate(this%m_work(m, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_work(m, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_ywork)) then
        if (size(this%m_ywork) /= neqn) then
            deallocate(this%m_ywork)
            allocate(this%m_ywork(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_ywork(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_dydx)) then
        if (size(this%m_dydx) /= neqn) then
            deallocate(this%m_dydx)
            allocate(this%m_dydx(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_dydx(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_jac)) then
        if (size(this%m_jac, 1) /= neqn .or. size(this%m_jac, 2) /= neqn) then
            deallocate(this%m_jac)
            allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_dfdx)) then
        if (size(this%m_dfdx) /= neqn) then
            deallocate(this%m_dfdx)
            allocate(this%m_dfdx(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_dfdx(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_a)) then
        if (size(this%m_a, 1) /= neqn .or. size(this%m_a, 2) /= neqn) then
            deallocate(this%m_a)
            allocate(this%m_a(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_a(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_iwork)) then
        if (size(this%m_iwork) /= neqn) then
            deallocate(this%m_iwork)
            allocate(this%m_iwork(neqn), stat = flag, source = 0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_iwork(neqn), stat = flag, source = 0)
        if (flag /= 0) go to 10
    end if
    
    if (allocated(this%m_interp)) then
        if (size(this%m_interp, 1) /= neqn .or. &
            size(this%m_interp, 2) /= 4) &
        then
            deallocate(this%m_interp)
            allocate(this%m_interp(neqn, 4), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_interp(neqn, 4), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("rbk_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, xmax
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: neqn
    real(real64) :: h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)

    ! Initialize the workspaces
    call this%allocate_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! If this is the first step, we need to compute a function evaluation; else,
    ! we can rely upon the output from the previous step
    if (this%m_firstStep) then
        call sys%ode(x, y, this%m_dydx)
    end if

    ! Compute the Jacobian
    call sys%compute_jacobian(x, y, this%m_dydx, this%m_jac, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute df/dx.  Use a forward difference estimate such that the previous
    ! function evaluation can be used
    h = sys%get_finite_difference_step()
    call sys%ode(x + h, y, this%m_dfdx)
    this%m_dfdx = (this%m_dfdx - this%m_dydx) / h

    ! Call the base routine to move forward with the remainder of the step
    call vsi_step(this, sys, x, xmax, y, yn, err = errmgr)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
    fprev, err)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn, en
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err
 
    ! Local Variables
    integer(int32) :: i, j, neqn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)

    ! Set up the system of linear equations and compute its LU factorization
    do j = 1, neqn
        do i = 1, neqn
            if (i == j) then
                this%m_a(i,j) = 1.0d0 / (h * gamma) - this%m_jac(i,j)
            else
                this%m_a(i,j) = -this%m_jac(i,j)
            end if
        end do
    end do
    call lu_factor(this%m_a, this%m_iwork, errmgr)
    if (errmgr%has_error_occurred()) return

    ! K1
    this%m_work(:,1) = this%m_dydx + h * d1 * this%m_dfdx
    call solve_lu(this%m_a, this%m_iwork, this%m_work(:,1))

    ! K2
    this%m_work(:,2) = y + a21 * this%m_work(:,1)
    call sys%ode(x + c2 * h, this%m_work(:,2), this%m_ywork)
    this%m_work(:,2) = this%m_ywork + h * d2 * this%m_dfdx + &
        c21 * this%m_work(:,1) / h
    call solve_lu(this%m_a, this%m_iwork, this%m_work(:,2))

    ! K3
    this%m_work(:,3) = y + a31 * this%m_work(:,1) + a32 * this%m_work(:,2)
    call sys%ode(x + c3 * h, this%m_work(:,3), this%m_ywork)
    this%m_work(:,3) = this%m_ywork + h * d3 * this%m_dfdx + &
        (c31 * this%m_work(:,1) + c32 * this%m_work(:,2)) / h
    call solve_lu(this%m_a, this%m_iwork, this%m_work(:,3))

    ! K4
    this%m_work(:,4) = y + a41 * this%m_work(:,1) + a42 * this%m_work(:,2) + &
        a43 * this%m_work(:,3)
    call sys%ode(x + c4 * h, this%m_work(:,4), this%m_ywork)
    this%m_work(:,4) = this%m_ywork + h * d4 * this%m_dfdx + &
        (c41 * this%m_work(:,1) + c42 * this%m_work(:,2) + &
        c43 * this%m_work(:,3)) / h
    call solve_lu(this%m_a, this%m_iwork, this%m_work(:,4))

    ! K5
    this%m_work(:,6) = y + a51 * this%m_work(:,1) + a52 * this%m_work(:,2) + &
        a53 * this%m_work(:,3) + a54 * this%m_work(:,4)
    call sys%ode(x + h, this%m_work(:,6), this%m_ywork)
    this%m_work(:,5) = this%m_ywork + (c51 * this%m_work(:,1) + &
        c52 * this%m_work(:,2) + c53 * this%m_work(:,3) + &
        c54 * this%m_work(:,4)) / h
    call solve_lu(this%m_a, this%m_iwork, this%m_work(:,5))

    ! Compute the error estimate
    this%m_work(:,6) = this%m_work(:,6) + this%m_work(:,5)
    call sys%ode(x + h, this%m_work(:,6), this%m_ywork)
    en = this%m_ywork + (c61 * this%m_work(:,1) + c62 * this%m_work(:,2) + &
        c3 * this%m_work(:,3) + c4 * this%m_work(:,4) + &
        c5 * this%m_work(:,5)) / h
    call solve_lu(this%m_a, this%m_iwork, en)   ! error estimate

    ! Compute the new solution estimate
    yn = this%m_work(:,6) + en

    ! End
    return
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_reset(this)
    class(rosenbrock_integrator), intent(inout) :: this
    this%m_firstStep = .true.
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_on_successful_step(this, x, xn, y, yn)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn

    print *, ""
    print *, x
    print *, this%get_step_size()
    print *, this%m_enormPrev

    ! Set up the interpolation polynomial
    call this%set_up_interpolation(x, xn, y, yn, this%m_work)

    ! Store the last results as the first
    this%m_firstStep = .false.
    this%m_dydx = this%m_ywork
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_set_up_interp(this, x, xn, y, yn, k)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn
    real(real64), intent(in), dimension(:,:) :: k

    ! Process
    this%m_interp(:,1) = y
    this%m_interp(:,2) = yn
    this%m_interp(:,3) = d21 * k(:,1) + d22 * k(:,2) + d23 * k(:,3) + &
        d24 * k(:,4) + d25 * k(:,5)
    this%m_interp(:,4) = d31 * k(:,1) + d32 * k(:,2) + d33 * k(:,3) + &
        d34 * k(:,4) + d35 * k(:,5)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rbk_interp(this, xprev, xnew, x, y, err)
    ! Arguments
    class(rosenbrock_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err
 
    ! Local Variables
    integer(int32) :: neqn
    real(real64) :: s, s1, h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(this%m_work, 1)
    h = xnew - xprev
    s = (x - xprev) / h
    s1 = 1.0d0 - s

    ! Input Check
    if (size(y) /= neqn) go to 10

    ! Process
    y = this%m_interp(:,1) * s + s * (this%m_interp(:,2) + &
        s1 * (this%m_interp(:,3) + s * this%m_interp(:,4)))

    ! End
    return
 
    ! Y is not sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to be of length ", &
        neqn, ", but was found to be of length ", size(y), "."
    call errmgr%report_error("bsrk32_interp", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return
 
    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
pure module function rbk_get_order(this) result(rst)
    class(rosenbrock_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
module function rbk_next_step(this, hn, en, enm1) result(rst)
    ! Arguments
    class(rosenbrock_integrator), intent(inout) :: this
    real(real64), intent(in) :: hn, en, enm1
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: fac1 = 5.0d0
    real(real64), parameter :: fac2 = 1.0d0 / 6.0d0
    
    ! Local Variables
    integer(int32) :: k
    real(real64) :: s, fac, facpred, hnew, maxstep, eold, hold

    ! Process
    k = this%get_order()
    s = this%get_safety_factor()
    fac = max(fac2, min(fac1, en**(1.0d0 / k)))
    eold = max(1.0d-2, enm1)
    hold = this%m_prevStep
    hnew = hn / fac
    if (en <= 1.0d0) then
        if (.not.this%m_firstStep) then
            facpred = (hold / hn) * (en**2 / eold)**(1.0d0 / k)
            facpred = max(fac2, min(fac1, facpred))
            fac = max(fac, facpred)
            hnew = hn / fac
        end if
        rst = s * hnew
        this%m_prevStep = hn
    else
        rst = s * hnew
    end if

    maxstep = abs(this%get_max_step_size())
    if (abs(rst) > maxstep) rst = sign(maxstep, rst)
end function

! ------------------------------------------------------------------------------
end submodule