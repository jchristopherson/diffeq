submodule (diffeq) diffeq_rkvariable
contains
! ------------------------------------------------------------------------------
module subroutine rkv_alloc_workspace(this, neqn, err)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: m, n, flag
    character(len = :), allocatable :: errmsg

    ! Process
    m = neqn
    n = this%get_stage_count()
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

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("rkv_alloc_workspace", trim(errmsg), &
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
    call this%define_model()

    ! Ensure the workspace arrays are allocated
    call this%allocate_rkv_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! The Butcher tableau is lower triangular as this is an explicit integrator
    if (.not.this%is_fsal() .or. this%m_firstStep) then
        ! On FSAL integrators, we only need to make this call on the first step
        ! as the integrator uses the last evaluation from the previous step
        ! as this step.  On non-FSAL integrators we always need to compute an
        ! updated first step.
        call sys%fcn(x, y, this%m_work(:,1))
    end if
    do i = 2, n
        this%m_ywork = 0.0d0
        do j = 1, i - 1 ! only reference the sub-diagonal components
            this%m_ywork = this%m_ywork + this%get_method_factor(i, j) * &
                this%m_work(:,j)
        end do

        call sys%fcn( &
            x + h * this%get_position_factor(i), &
            y + h * this%m_ywork, &
            this%m_work(:,i) &  ! output
        )
    end do

    ! Compute the two solution estimates, and the resulting error estimate
    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%get_quadrature_weight(i) * this%m_work(:,i)
        else
            this%m_ywork = this%m_ywork + this%get_quadrature_weight(i) * &
                this%m_work(:,i)
        end if
    end do
    yn = y + h * this%m_ywork

    do i = 1, n
        if (i == 1) then
            this%m_ywork = this%get_error_factor(i) * this%m_work(:,i)
        else
            this%m_ywork = this%m_ywork + this%get_error_factor(i) * &
                this%m_work(:,i)
        end if
    end do
    en = h * this%m_ywork
end subroutine

! ------------------------------------------------------------------------------
module subroutine rkv_on_successful_step(this)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this

    ! Local Variables
    integer(int32) :: n

    ! Store the last result as the first, if this is FSAL
    if (this%is_fsal()) then
        this%m_firstStep = .false.
        n = this%get_stage_count()
        this%m_work(:,1) = this%m_work(:,n)
    end if
end subroutine

! ------------------------------------------------------------------------------
end submodule