! This file is part of diffeq.
! 
! diffeq is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! diffeq is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with diffeq. If not, see <https://www.gnu.org/licenses/>.

program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    implicit none

    ! Initial Conditions & Time Constraints
    real(real64), parameter :: t(2) = [0.0d0, 5.0d1]
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]

    ! Local Variables
    type(runge_kutta_23) :: integrator_1
    type(runge_kutta_45) :: integrator_2
    type(runge_kutta_853) :: integrator_3
    type(rosenbrock) :: integrator_4
    type(kennedy_carpenter_4) :: integrator_5
    type(kennedy_carpenter_5) :: integrator_6
    type(tsitouras_54) :: integrator_7
    type(bdf) :: integrator_8
    type(adams) :: integrator_9
    type(ode_container) :: mdl
    type(ode_results) :: s1, s2, s2a, s3, s4, s5, s6, s7, s8, s9

    ! Define the model
    mdl%fcn => vanderpol

    ! Integrate the model with each integrator
    call integrator_1%solve(mdl, t, ic)
    call integrator_2%solve(mdl, t, ic)
    call integrator_3%solve(mdl, t, ic)
    call integrator_4%solve(mdl, t, ic)
    call integrator_5%solve(mdl, t, ic)
    call integrator_6%solve(mdl, t, ic)
    call integrator_7%solve(mdl, t, ic)
    call integrator_8%solve(mdl, t, ic)
    call integrator_9%solve(mdl, t, ic)

    ! Retrieve the solution from each integrator
    s1 = integrator_1%get_ode_results()
    s2 = integrator_2%get_ode_results()
    s3 = integrator_3%get_ode_results()
    s4 = integrator_4%get_ode_results()
    s5 = integrator_5%get_ode_results()
    s6 = integrator_6%get_ode_results()
    s7 = integrator_7%get_ode_results()
    s8 = integrator_8%get_ode_results()
    s9 = integrator_9%get_ode_results()

    ! Print out the size of each solution
    print "(A)", "RUNGE_KUTTA_23: " 
    print "(A, I0, A)", achar(9), size(s1, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s1%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s1%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s1%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "RUNGE_KUTTA_45: "
    print "(A, I0, A)", achar(9), size(s2, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s2%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s2%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s2%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "RUNGE_KUTTA_853: "
    print "(A, I0, A)", achar(9), size(s3, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s3%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s3%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s3%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "ROSENBROCK: "
    print "(A, I0, A)", achar(9), size(s4, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s4%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s4%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s4%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s4%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s4%statistics%matrix_factorizations, " Matrix Factorizations"

    ! Now, implement a PI controller and check its effect.  This will likely
    ! increase the number of steps (loss of efficiency), but if there were
    ! any stability issues, stability will likely improve.  Stability is likely
    ! not relevant on this problem, but it's here for illustration purposes.
    call integrator_2%clear_buffer()
    call integrator_2%set_step_size_control_parameter(0.1d0)
    call integrator_2%solve(mdl, t, ic)
    s2a = integrator_2%get_ode_results()
    print "(A)", "RUNGE_KUTTA_45 w/ PI Controller: "
    print "(A, I0, A)", achar(9), size(s2a, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s2a%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s2a%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s2a%statistics%function_evaluations, " Function Evaluations"

    ! Kennedy-Carpenter Integrators
    print "(A)", "KC4: "
    print "(A, I0, A)", achar(9), size(s5, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s5%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s5%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s5%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s5%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s5%statistics%matrix_factorizations, " Matrix Factorizations"

    print "(A)", "KC5: "
    print "(A, I0, A)", achar(9), size(s6, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s6%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s6%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s6%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s6%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s6%statistics%matrix_factorizations, " Matrix Factorizations"

    ! Additional Integrators
    print "(A)", "TSITOURAS 4/5: "
    print "(A, I0, A)", achar(9), size(s7, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s7%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s7%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s7%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "BDF: "
    print "(A, I0, A)", achar(9), size(s8, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s8%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s8%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s8%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s8%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s8%statistics%matrix_factorizations, " Matrix Factorizations"

    print "(A)", "ADAMS: "
    print "(A, I0, A)", achar(9), size(s9, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s9%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s9%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s9%statistics%function_evaluations, " Function Evaluations"
end program