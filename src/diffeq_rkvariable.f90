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
module subroutine rkv_attempt_step(this, sys, h, x, y, yn, en)
    ! Arguments
    class(rk_variable_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn, en

    ! Local Variables
    integer(int32) :: i, j, n

    ! Initialization
    n = this%get_stage_count()

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
end subroutine

! ------------------------------------------------------------------------------
end submodule