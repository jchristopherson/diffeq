submodule (diffeq) diffeq_implicit_rk
contains
! ------------------------------------------------------------------------------
pure module function irk_get_old_step(this) result(rst)
    class(implicit_rk_variable_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_hold
end function

! ---------------------
module subroutine irk_set_old_step(this, x)
    class(implicit_rk_variable_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_hold = x
end subroutine

! ------------------------------------------------------------------------------
module function irk_compute_next_step_size(this, hn, en, enm1) result(rst)
    ! Arguments
    class(implicit_rk_variable_integrator), intent(inout) :: this
    real(real64), intent(in) :: hn, en, enm1
    real(real64) :: rst

    ! Parameters
    real(real64), parameter :: fac1 = 5.0d0
    real(real64), parameter :: fac2 = 1.0d0 / 8.0d0

    ! Local Variables
    real(real64) :: fac, safe, hnew, facpred, e, hold, hmax

    ! Initialization
    safe = this%get_safety_factor()
    e = 1.0d0 / real(this%get_order(), real64)
    fac = max(fac2, min(fac1, (en**2 / enm1)**e / safe))
    rst = hn / fac
    hold = this%get_previous_step_size()
    hmax = this%get_max_step_size()

    ! Process
    if (en <= 1.0d0) then
        if (.not.this%m_firstStep) then
            facpred = (hold / hn) * (err**2 / enm1)**e / safe
            facpred = max(fac2, min(fac1, facpred))
            fac = max(fac, facpred)
            rst = hn / fac
        end if
        this%m_firstStep = .false.
        call this%set_previous_step_size(hn)
        this%m_enormPrev = max(1.0d-2, en)
    end if
    if (abs(rst) > abs(hmax)) rst = sign(hmax, rst)
end function

! ------------------------------------------------------------------------------
pure module function irk_get_is_jac_current(this) result(rst)
    class(implicit_rk_variable_integrator), intent(in) :: this
    logical :: rst
    rst = this%m_isJacCurrent
end function

! ---------------------
module subroutine irk_set_is_jac_current(this, x)
    class(implicit_rk_variable_integrator), intent(inout) :: this
    logical, intent(in) :: x
    this%m_isJacCurrent = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine irk_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(implicit_rk_variable_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, xmax
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    real(real64) :: h
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    h = this%get_next_step_size()

    ! Compute the system matrices
    call this%build_factored_newton_matrix(sys, h, x, y, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Call the base routine to finish the step
    call vsi_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, errmgr)
end subroutine

! ------------------------------------------------------------------------------
pure module function irk_get_max_newton_iter(this) result(rst)
    class(implicit_rk_variable_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_maxNewtonIter
end function

! --------------------
module subroutine irk_set_max_newton_iter(this, x)
    class(implicit_rk_variable_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_maxNewtonIter = x
end subroutine

! ------------------------------------------------------------------------------
end submodule