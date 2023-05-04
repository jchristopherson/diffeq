submodule (diffeq) diffeq_dprk45
    ! Dormand-Prince 4th/5th Order Model Coefficients
    real(real64), parameter :: a21 = 1.0d0 / 5.0d0
    real(real64), parameter :: a31 = 3.0d0 / 40.0d0
    real(real64), parameter :: a32 = 9.0d0 / 40.0d0
    real(real64), parameter :: a41 = 44.0d0 / 45.0d0
    real(real64), parameter :: a42 = -56.0d0 / 15.0d0
    real(real64), parameter :: a43 = 32.0d0 / 9.0d0
    real(real64), parameter :: a51 = 1.9372d4 / 6.561d3
    real(real64), parameter :: a52 = -2.536d4 / 2.187d3
    real(real64), parameter :: a53 = 6.4448d4 / 6.561d3
    real(real64), parameter :: a54 = -2.12d2 / 7.29d2
    real(real64), parameter :: a61 = 9.017d3 / 3.168d3
    real(real64), parameter :: a62 = -3.55d2 / 33.0d0
    real(real64), parameter :: a63 = 4.6732d4 / 5.247d3
    real(real64), parameter :: a64 = 49.0d0 / 1.76d2
    real(real64), parameter :: a65 = -5.103d3 / 1.8656d4
    real(real64), parameter :: a71 = 35.0d0 / 3.84d2
    real(real64), parameter :: a72 = 0.0d0
    real(real64), parameter :: a73 = 5.0d2 / 1.113d3
    real(real64), parameter :: a74 = 1.25d2 / 1.92d2
    real(real64), parameter :: a75 = -2.187d3 / 6.784d3
    real(real64), parameter :: a76 = 11.0d0 / 84.0d0
    
    real(real64), parameter :: e1 = -71.0d0 / 5.76d4
    real(real64), parameter :: e2 = 0.0d0
    real(real64), parameter :: e3 = 71.0d0 / 1.6695d4
    real(real64), parameter :: e4 = -71.0d0 / 1.92d3
    real(real64), parameter :: e5 = 1.7253d4 / 3.392d5
    real(real64), parameter :: e6 = -22.0d0 / 5.25d2
    real(real64), parameter :: e7 = 1.0d0 / 4.0d1

    real(real64), parameter :: c2 = 1.0d0 / 5.0d0
    real(real64), parameter :: c3 = 3.0d0 / 1.0d1
    real(real64), parameter :: c4 = 4.0d0 / 5.0d0
    real(real64), parameter :: c5 = 8.0d0 / 9.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: c7 = 1.0d0
contains
! ------------------------------------------------------------------------------
module subroutine dprk45_define_model(this)
    ! Arguments
    class(dprk45_integrator), intent(inout) :: this

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

    this%m_a(5,1) = a51
    this%m_a(5,2) = a52
    this%m_a(5,3) = a53
    this%m_a(5,4) = a54

    this%m_a(6,1) = a61
    this%m_a(6,2) = a62
    this%m_a(6,3) = a63
    this%m_a(6,4) = a64
    this%m_a(6,5) = a65
    
    this%m_a(7,1) = a71
    this%m_a(7,2) = a72
    this%m_a(7,3) = a73
    this%m_a(7,4) = a74
    this%m_a(7,5) = a75
    this%m_a(7,6) = a76

    ! B
    this%m_b(1) = a71
    this%m_b(2) = a72
    this%m_b(3) = a73
    this%m_b(4) = a74
    this%m_b(5) = a75
    this%m_b(6) = a76
    this%m_b(7) = 0.0d0

    ! C
    this%m_c(1) = 0.0d0
    this%m_c(2) = c2
    this%m_c(3) = c3
    this%m_c(4) = c4
    this%m_c(5) = c5
    this%m_c(6) = c6
    this%m_c(7) = c7

    ! E
    this%m_e(1) = e1
    this%m_e(2) = e2
    this%m_e(3) = e3
    this%m_e(4) = e4
    this%m_e(5) = e5
    this%m_e(6) = e6
    this%m_e(7) = e7

    ! Update definition status
    this%m_modelDefined = .true.
end subroutine

! ------------------------------------------------------------------------------
pure module function dprk45_get_method_factor(this, i, j) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32), intent(in) :: i, j
    real(real64) :: rst
    rst = this%m_a(i, j)
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_quad_weights(this, i) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_b(i)
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_error_factor(this, i) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_e(i)
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_position_factor(this, i) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    rst = this%m_c(i)
end function

! ------------------------------------------------------------------------------
pure module function dprk45_is_fsal(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_stage_count(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 7
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_order(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 5
end function

! ------------------------------------------------------------------------------
end submodule