submodule (diffeq) diffeq_sdirk4
    ! Model Parameters (ESDIRK4(3)6L[2]SA, Table 16, pg 90)
    ! https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf
    real(real64), parameter :: gamma = 0.25d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = (2.0d0 - sqrt(2.0d0)) / 4.0d0
    real(real64), parameter :: c4 = 5.0d0 / 8.0d0
    real(real64), parameter :: c5 = 26.0d0 / 25.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: a21 = gamma
    real(real64), parameter :: a22 = gamma
    real(real64), parameter :: a32 = (1.0d0 - sqrt(2.0d0)) / 8.0d0
    real(real64), parameter :: a31 = (c3 - a32 - gamma)
    real(real64), parameter :: a33 = gamma
    real(real64), parameter :: a43 = 7.0d0 * (1.0d0 + sqrt(2.0d0)) / 32.0d0
    real(real64), parameter :: a42 = (5.0d0 - 7.0d0 * sqrt(2.0d0)) / 64.0d0
    real(real64), parameter :: a41 = c4 - a42 - a43 - gamma
    real(real64), parameter :: a44 = gamma
    real(real64), parameter :: a54 = &
        166.0d0 * (-97.0d0 + 376.0d0 * sqrt(2.0d0)) / 1.09375d5
    real(real64), parameter :: a53 = &
        (5.06605d5 + 1.32109d5 * sqrt(2.0d0)) / 4.375d5
    real(real64), parameter :: a52 = &
        (-1.3796d4 - 5.4539d4 * sqrt(2.0d0)) / 1.25d5
    real(real64), parameter :: a51 = c5 - a52 - a53 - a54 - gamma
    real(real64), parameter :: a55 = gamma
    real(real64), parameter :: b5 = &
        -15.625d3 * (97.0d0 + 376.0d0 * sqrt(2.0d0)) / 9.0749876d7
    real(real64), parameter :: b4 = &
        -16.0d0 * (-2.2922d4 + 3.525d3 * sqrt(2.0d0)) / 5.71953d5
    real(real64), parameter :: b3 = &
        47.0d0 * (-267.0d0 + 1.783d3 * sqrt(2.0d0)) / 2.73343d5
    real(real64), parameter :: b2 = (1.181d3 - 9.87d2 * sqrt(2.0d0)) / 1.3782d4
    real(real64), parameter :: b1 = 1.0d0 - b2 - b3 - b4 - b5 - gamma
    real(real64), parameter :: b6 = gamma
    real(real64), parameter :: b2a = -4.80923228411d11 / 4.982971448372d12
    real(real64), parameter :: b3a = 6.709447293961d12 / 1.2833189095359d13
    real(real64), parameter :: b4a = 3.513175791894d12 / 6.748737351361d12
    real(real64), parameter :: b5a = -4.98863281070d11 / 6.042575550617d12
    real(real64), parameter :: b6a = 2.077005547802d12 / 8.945017530137d12
    real(real64), parameter :: b1a = 1.0d0 - b2a - b3a - b4a - b5a - b6a

contains
! ------------------------------------------------------------------------------
pure module function sd4_get_order(this) result(rst)
    class(sdirk4_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure module function sd4_get_stage_count(this) result(rst)
    class(sdirk4_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 6
end function

! ------------------------------------------------------------------------------
module subroutine sd4_define_model(this)
    ! Arguments
    class(sdirk4_integrator), intent(inout) :: this

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%a = 0.0d0
    
    this%a(1,1) = a22

    this%a(2,1) = a21
    this%a(2,2) = a22

    this%a(3,1) = a31
    this%a(3,2) = a32
    this%a(3,3) = a33

    this%a(4,1) = a41
    this%a(4,2) = a42
    this%a(4,3) = a43
    this%a(4,4) = a44

    this%a(5,1) = a51
    this%a(5,2) = a52
    this%a(5,3) = a53
    this%a(5,4) = a54
    this%a(5,5) = a55

    this%a(6,1) = b1
    this%a(6,2) = b2
    this%a(6,3) = b3
    this%a(6,4) = b4
    this%a(6,5) = b5
    this%a(6,6) = b6

    ! B
    this%b(1) = b1
    this%b(2) = b2
    this%b(3) = b3
    this%b(4) = b4
    this%b(5) = b5
    this%b(6) = b6

    ! C
    this%c(1) = 0.0d0
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4
    this%c(5) = c5
    this%c(6) = c6

    ! E
    this%e(1) = b1a - b1
    this%e(2) = b2a - b2
    this%e(3) = b3a - b3
    this%e(4) = b4a - b4
    this%e(5) = b5a - b5
    this%e(6) = b6a - b6

    ! Update definition status
    this%m_modelDefined = .true.
end subroutine

! ------------------------------------------------------------------------------
pure module function sd4_is_fsal(this) result(rst)
    class(sdirk4_integrator), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
module subroutine sd4_set_up_interp(this, x, xn, y, yn, k)
    ! Arguments
    class(sdirk4_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn
    real(real64), intent(in), dimension(:,:) :: k

    ! Local Variables
end subroutine

! ------------------------------------------------------------------------------
module subroutine sd4_interp(this, xprev, xnew, x, y, err)
    ! Arguments
    class(sdirk4_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
end subroutine

! ------------------------------------------------------------------------------
end submodule