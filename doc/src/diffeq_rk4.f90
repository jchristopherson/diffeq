module diffeq_rk4
    !! This module provides a 4th order Runge-Kutta fixed-step integrator.
    use iso_fortran_env
    use diffeq_rk_fixed_step
    implicit none
    private
    public :: rk4_fixed_integrator

    type, extends(rk_fixed_integrator) :: rk4_fixed_integrator
        !! Defines the explicit, 4th order, Runge-Kutta fixed-step 
        !! integrator.
    contains
        procedure, public :: get_order => rk4_get_order
            !! Returns the order of the integrator.
        procedure, public :: get_method_factor => rk4_get_a_param
            !! Gets the requested method factor from the Butcher tableau.
        procedure, public :: get_quadrature_weight => rk4_get_b_param
            !! Gets the requested quadrature weight from the Butcher tableau.
        procedure, public :: get_position_factor => rk4_get_c_param
            !! Gets the requested position factor from the Butcher tableau.
    end type

contains
! ------------------------------------------------------------------------------
pure function rk4_get_order(this) result(rst)
    !! Returns the order of the integrator.
    class(rk4_fixed_integrator), intent(in) :: this
        !! The rk4_fixed_integrator object.
    integer(int32) :: rst
        !! The order of the integrator.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function rk4_get_a_param(this, i, j) result(rst)
    !! Gets the requested method factor from the Butcher tableau.
    class(rk4_fixed_integrator), intent(in) :: this
        !! The rk4_fixed_integrator object.
    integer(int32), intent(in) :: i
        !! The row index of the parameter from the Butcher tableau.
    integer(int32), intent(in) :: j
        !! The column index of the parameter from the Butcher tableau.
    real(real64) :: rst
        !! The requested parameter.
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
pure function rk4_get_b_param(this, i) result(rst)
    !! Gets the requested quadrature weight from the Butcher tableau.
    class(rk4_fixed_integrator), intent(in) :: this
        !! The rk4_fixed_integrator object.
    integer(int32), intent(in) :: i
        !! The index of the parameter from the Butcher tableau.
    real(real64) :: rst
        !! The requested parameter.
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
pure function rk4_get_c_param(this, i) result(rst)
    !! Gets the requested position factor from the Butcher tableau.
    class(rk4_fixed_integrator), intent(in) :: this
        !! The rk4_fixed_integrator object.
    integer(int32), intent(in) :: i
        !! The index of the parameter from the Butcher tableau.
    real(real64) :: rst
        !! The requested parameter.
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
end module