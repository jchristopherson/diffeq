submodule (diffeq) diffeq_rk4
contains
! ------------------------------------------------------------------------------
pure module function rk4_get_order(this) result(rst)
    class(rk4_fixed_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure module function rk4_get_a_param(this, i, j) result(rst)
    class(rk4_fixed_integrator), intent(in) :: this
    integer(int32), intent(in) :: i, j
    real(real64) :: rst
    if (i == 2 .and. j == 1) then
        rst = 0.5d0
    else if (i == 3 .and. j == 2) then
        rst = 0.5d0
    else if (i == 4 .and. j == 3) then
        rst = 1.0d0
    else
        rst = 0.0d0
    end if
end function

! ------------------------------------------------------------------------------
pure module function rk4_get_b_param(this, i) result(rst)
    class(rk4_fixed_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    real(real64), parameter :: sixth = 1.0d0 / 6.0d0
    real(real64), parameter :: third = 1.0d0 / 3.0d0
    select case (i)
    case (1)
        rst = sixth
    case (2)
        rst = third
    case (3)
        rst = third
    case (4)
        rst = sixth
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
pure module function rk4_get_c_param(this, i) result(rst)
    class(rk4_fixed_integrator), intent(in) :: this
    integer(int32), intent(in) :: i
    real(real64) :: rst
    select case (i)
    case (2)
        rst = 0.5d0
    case (3)
        rst = 0.5d0
    case default
        rst = 0.0d0
    end select
end function

! ------------------------------------------------------------------------------
end submodule