submodule (diffeq) diffeq_vs_integrator
contains
! ------------------------------------------------------------------------------
pure module function vsi_get_safety_factor(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_safetyfactor
end function

! --------------------
module subroutine vsi_set_safety_factor(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_safetyfactor = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_alpha(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_alpha
end function

! --------------------
module subroutine vsi_set_alpha(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64) :: x
    this%m_alpha = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_beta(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_beta
end function

! --------------------
module subroutine vsi_set_beta(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_beta = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_get_max_step(this) result(rst)
    class(variable_step_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_maxstep
end function

! --------------------
module subroutine vsi_set_max_step(this, x)
    class(variable_step_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_maxstep = x
end subroutine

! ------------------------------------------------------------------------------
pure module function vsi_estimate_error(this, y, ys, atol, rtol) result(rst)
    ! Arguments
    class(variable_step_integrator), intent(in) :: this
    real(real64), intent(in), dimension(:) :: y, ys, atol, rtol
    real(real64) :: rst
    
    ! Process
    rst = norm2((y - ys) / (atol + maxval(abs(y)) * rtol))
end function

! ------------------------------------------------------------------------------
pure module function vsi_next_step(this, hn, en, enm1) result(rst)
    ! Arguments
    class(variable_step_integrator), intent(in) :: this
    real(real64), intent(in) :: hn, en, enm1
    real(real64) :: rst

    ! Arguments
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
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule