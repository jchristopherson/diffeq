submodule (diffeq) diffeq_rkvariable
contains
! ------------------------------------------------------------------------------
module subroutine rkv_alloc_workspace(this, neqn, err)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: m, n, flag
    character(len = :), allocatable :: errmsg
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Process
    m = neqn
    n = this%get_stage_count()
    if (allocated(this%f)) then
        if (size(this%f, 1) /= m .or. size(this%f, 2) /= n) then
            deallocate(this%f)
            allocate(this%f(m, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%f(m, n), stat = flag, source = 0.0d0)
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

    if (allocated(this%a)) then
        if (size(this%a, 1) /= n .or. size(this%a, 2) /= n) then
            deallocate(this%a)
            allocate(this%a(n, n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%a(n, n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%b)) then
        if (size(this%b) /= n) then
            deallocate(this%b)
            allocate(this%b(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%b(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%c)) then
        if (size(this%c) /= n) then
            deallocate(this%c)
            allocate(this%c(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%c(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%e)) then
        if (size(this%e) /= n) then
            deallocate(this%e)
            allocate(this%e(n), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%e(n), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! Define the model parameters
    call this%define_model()

    ! Call the base method
    call vsi_alloc_workspace(this, neqn, errmgr)

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("rkv_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rkv_reset(this)
    class(rk_variable_integrator), intent(inout) :: this
    this%m_firstStep = .true.
end subroutine

! ------------------------------------------------------------------------------
module subroutine rkv_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
    fprev, err)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn, en
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, j, n, neqn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = this%get_stage_count()
    neqn = size(y)

    ! The Butcher tableau is lower triangular as this is an explicit integrator
    if (.not.this%is_fsal() .or. this%m_firstStep) then
        ! On FSAL integrators, we only need to make this call on the first step
        ! as the integrator uses the last evaluation from the previous step
        ! as this step.  On non-FSAL integrators we always need to compute an
        ! updated first step.
        call sys%ode(x, y, this%f(:,1))
    end if
    do i = 2, n
        this%m_ywork = 0.0d0
        do j = 1, i - 1 ! only reference the sub-diagonal components
            this%m_ywork = this%m_ywork + this%a(i, j) * &
                this%f(:,j)
        end do

        call sys%ode( &
            x + h * this%c(i), &
            y + h * this%m_ywork, &
            this%f(:,i) &  ! output
        )
    end do

    ! Compute the two solution estimates, and the resulting error estimate
    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%b(i) * this%f(:,i)
        else
            this%m_ywork = this%m_ywork + this%b(i) * &
                this%f(:,i)
        end if
    end do
    yn = y + h * this%m_ywork

    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%e(i) * this%f(:,i)
        else
            this%m_ywork = this%m_ywork + this%e(i) * &
                this%f(:,i)
        end if
    end do
    en = h * this%m_ywork
end subroutine

! ------------------------------------------------------------------------------
module subroutine rkv_on_successful_step(this, x, xn, y, yn)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn

    ! Local Variables
    integer(int32) :: n

    ! Set up the interpolation polynomial - TO DO check for dense output first
    call this%set_up_interpolation(x, xn, y, yn, this%f)

    ! Store the last result as the first, if this is FSAL
    if (this%is_fsal()) then
        this%m_firstStep = .false.
        n = this%get_stage_count()
        this%f(:,1) = this%f(:,n)
    end if
end subroutine


! ------------------------------------------------------------------------------
pure module function rkv_get_alpha(this) result(rst)
    class(rk_variable_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_alpha
end function

! --------------------
module subroutine rkv_set_alpha(this, x)
    class(rk_variable_integrator), intent(inout) :: this
    real(real64) :: x
    this%m_alpha = x
end subroutine

! ------------------------------------------------------------------------------
pure module function rkv_get_beta(this) result(rst)
    class(rk_variable_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_beta
end function

! --------------------
module subroutine rkv_set_beta(this, x)
    class(rk_variable_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_beta = x
end subroutine

! ------------------------------------------------------------------------------
module function rkv_next_step(this, hn, en, enm1) result(rst)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    real(real64), intent(in) :: hn, en, enm1
    real(real64) :: rst

    ! Local Variables
    integer(int32) :: k
    real(real64) :: s, hest, a, b, maxstep

    ! Process
    k = this%get_order()
    s = this%get_safety_factor()
    hest = s * hn * (1.0d0 / en)**(1.0d0 / k)
    
    a = this%get_alpha() / k
    b = this%get_beta() / k
    rst = hest * (en**a) * (enm1**b)

    maxstep = abs(this%get_max_step_size())
    if (abs(rst) > maxstep) rst = sign(maxstep, rst)
    if (rst / hn > 2.0d0) rst = 2.0d0 * hn
end function

! ------------------------------------------------------------------------------
end submodule