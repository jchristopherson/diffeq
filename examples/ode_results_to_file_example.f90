program example
    use iso_fortran_env
    use diffeq
    implicit none

    type(ode_container) :: model
    type(tsitouras_54) :: integrator
    type(ode_results) :: rst

    model%fcn => rhs
    call integrator%solve(model, [0.0d0, 0.5d0], [1.0d0, 0.0d0])
    rst = integrator%get_ode_results()
    call write_ode_results("ode_results_output.txt", rst)

contains
    subroutine rhs(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args

        real(real64), parameter :: zeta = 0.1d0
        real(real64), parameter :: wn = 1.0d2

        dxdt(1) = x(2)
        dxdt(2) = -2.0d0 * zeta * wn * x(2) - wn**2 * x(1)
    end subroutine
end program