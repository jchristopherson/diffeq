submodule (diffeq) diffeq_bsrk32
    real(real64), parameter :: a21 = 0.5d0
    real(real64), parameter :: a31 = 0.0d0
    real(real64), parameter :: a32 = 0.75d0
    real(real64), parameter :: a41 = 2.0d0 / 9.0d0
    real(real64), parameter :: a42 = 1.0d0 / 3.0d0
    real(real64), parameter :: a43 = 4.0d0 / 9.0d0

    real(real64), parameter :: b1 = 2.0d0 / 9.0d0
    real(real64), parameter :: b2 = 1.0d0 / 3.0d0
    real(real64), parameter :: b3 = 4.0d0 / 9.0d0
    real(real64), parameter :: b4 = 0.0d0

    real(real64), parameter :: b1a = 7.0d0 / 24.0d0
    real(real64), parameter :: b2a = 1.0d0 / 4.0d0
    real(real64), parameter :: b3a = 1.0d0 / 3.0d0
    real(real64), parameter :: b4a = 1.0d0 / 8.0d0

    real(real64), parameter :: c1 = 0.0d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = 0.75d0
    real(real64), parameter :: c4 = 1.0d0
contains
! ------------------------------------------------------------------------------
module subroutine bsrk32_define_model(this)
    ! Arguments
    class(bsrk32_integrator), intent(inout) :: this

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%m_a = 0.0d0

    this%m_a(2,1) = a21

    this%m_a(3,1) = a31
    this%m_a(3,2) = a32

    this%m_a(4,1) = a41
    this%m_a(4,2) = a42
    this%m_a(4,3) = a43

    ! B
    this%m_b(1) = b1
    this%m_b(2) = b2
    this%m_b(3) = b3
    this%m_b(4) = b4

    ! C
    this%m_c(1) = c1
    this%m_c(2) = c2
    this%m_c(3) = c3
    this%m_c(4) = c4

    ! E
    this%m_e(1) = b1a - b1
    this%m_e(2) = b2a - b2
    this%m_e(3) = b3a - b3
    this%m_e(4) = b4a - b4
end subroutine

! ------------------------------------------------------------------------------
pure module function bsrk32_get_method_factor(this, i, j) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32), intent(in) :: i, j
    real(real64) :: rst
    rst = this%m_a(i, j)
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_quad_weights(this, i) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_b(i)
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_error_factor(this, i) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_e(i)
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_position_factor(this, i) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_c(i)
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_is_fsal(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_stage_count(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_order(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 3
end function

! ------------------------------------------------------------------------------
module subroutine bsrk32_set_up_interp(this, y, yn, k)
    ! Arguments
    class(bsrk32_integrator), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: y, yn
    real(real64), intent(in), dimension(:,:) :: k
end subroutine

! ------------------------------------------------------------------------------
module subroutine bsrk32_interp(this, xprev, xnew, x, y, err)
    ! Arguments
    class(bsrk32_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule