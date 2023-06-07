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

! ------------------------------------------------------------------------------
end submodule