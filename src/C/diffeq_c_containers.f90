! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
module diffeq_c_containers
    use iso_c_binding
    use diffeq_containers
    use diffeq_container_io
    implicit none
    private
    public :: assignment(=)
    public :: c_ode_results_statistics

    interface assignment(=)
        module procedure :: c_ode_results_stats_assign
        module procedure :: ode_results_stats_to_c_assign
    end interface

    type, bind(C) :: c_ode_results_statistics
        !! A container of performance statistics for an
        !! integrator.
        integer(c_int) :: successful_steps = 0
            !! The number of successful steps taken by the integrator.
        integer(c_int) :: failed_steps = 0
            !! The number of unsuccessful steps taken by the integrator.
        integer(c_int) :: jacobian_evaluations = 0
            !! The number of Jacobian matrix evaluations.
        integer(c_int) :: function_evaluations = 0
            !! The number of function evaluations.
        integer(c_int) :: matrix_factorizations = 0
            !! The number of matrix factorizations.
    end type

    type, bind(C) :: c_ode_results
        !! A container for ODE results and integrator performance
        !! statistics.
        type(c_ptr) :: container
            !! A pointer to the underlying [[ode_results]] type.
        type(c_ode_results_statistics) :: statistics
            !! The integrator performance statistics.
    end type

contains
    subroutine c_ode_results_stats_assign(lhs, rhs)
        type(c_ode_results_statistics), intent(out) :: lhs
        type(ode_results_statistics) :: rhs
        lhs%successful_steps = integer(rhs%successful_steps, c_int)
        lhs%failed_steps = integer(rhs%failed_steps, c_int)
        lhs%jacobian_evaluations = integer(rhs%jacobian_evaluations, c_int)
        lhs%function_evaluations = integer(rhs%function_evaluations, c_int)
        lhs%matrix_factorizations = integer(rhs%matrix_factorizations, c_int)
    end subroutine

    subroutine ode_results_stats_to_c_assign(lhs, rhs)
        type(ode_results_statistics), intent(out) :: lhs
        type(c_ode_results_statistics), intent(in) :: rhs
        lhs%successful_steps = rhs%successful_steps
        lhs%failed_steps = rhs%failed_steps
        lhs%jacobian_evaluations = rhs%jacobian_evaluations
        lhs%function_evaluations = rhs%function_evaluations
        lhs%matrix_factorizations = rhs%matrix_factorizations
    end subroutine
end module