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
    use diffeq_c_errors
    use diffeq_c_helpers
    implicit none
    private
    public :: assignment(=)
    public :: c_ode_results_statistics
    public :: c_create_ode_results
    public :: c_free_ode_results

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

contains
! ------------------------------------------------------------------------------
    subroutine c_ode_results_stats_assign(lhs, rhs)
        type(c_ode_results_statistics), intent(out) :: lhs
        type(ode_results_statistics) :: rhs
        lhs%successful_steps = int(rhs%successful_steps, c_int)
        lhs%failed_steps = int(rhs%failed_steps, c_int)
        lhs%jacobian_evaluations = int(rhs%jacobian_evaluations, c_int)
        lhs%function_evaluations = int(rhs%function_evaluations, c_int)
        lhs%matrix_factorizations = int(rhs%matrix_factorizations, c_int)
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

! ------------------------------------------------------------------------------
    subroutine c_create_ode_results(c_obj, src)
        !! Constructs a new [[ode_results]] pointer and returns its address to
        !! the supplied C pointer.  Be sure to call c_free_ode_results once
        !! finished with this object.
        type(c_ptr), intent(out) :: c_obj
            !! A pointer to the newly created [[ode_results]] object.
        type(ode_results), intent(in) :: src
            !! The source [[ode_results]] object.

        type(ode_results), pointer :: obj
        allocate(obj, source = src)
        c_obj = c_loc(obj)
    end subroutine

    subroutine c_free_ode_results(c_obj) bind(C, name = "c_free_ode_results")
        !! Frees memory allocated by c_create_ode_results.
        type(c_ptr), value, intent(in) :: c_obj
            !! A pointer to the [[ode_results]] object to free.
        type(ode_results), pointer :: obj
        if (.not.c_associated(c_obj)) return
        call c_f_pointer(c_obj, obj)
        if (associated(obj)) deallocate(obj)
    end subroutine

    function c_get_ode_results_statistics(c_obj, stats) result(rst) &
        bind(C, name = "c_get_ode_results_statistics")
        !! Gets the integrator performance statistics from an [[ode_results]]
        !! object.
        type(c_ptr), value, intent(in) :: c_obj
            !! The pointer to the [[ode_results]] object.
        type(c_ode_results_statistics), intent(out) :: stats
            !! A [[c_ode_results_statistics]] container for the results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.  

        type(ode_results), pointer :: obj

        rst = C_NO_ERROR
        if (.not.c_associated(c_obj)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_pointer(c_obj, obj)
        stats = obj%statistics
    end function

    function c_get_ode_results_size(c_obj, rows, cols) result(rst) &
        bind(C, name = "c_get_ode_results_size")
        !! Gets the size of the matrix required to hold the integration results.
        type(c_ptr), value, intent(in) :: c_obj
            !! The pointer to the [[ode_results]] object.
        integer(c_int), intent(out) :: rows
            !! The number of rows in the matrix.
        integer(c_int), intent(out) :: cols
            !! The number of columns in the matrix.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        type(ode_results), pointer :: obj
        
        rst = C_NO_ERROR
        rows = 0
        cols = 0
        if (.not.c_associated(c_obj)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_pointer(c_obj, obj)
        rows = obj%count()
        cols = obj%state_variable_count() + 1
    end function

    function c_get_ode_results(c_obj, m, n, x, ldx) result(rst) &
        bind(C, name = "c_get_ode_results")
        !! Returns the contents of the [[ode_results]] object as a matrix.
        type(c_ptr), value, intent(in) :: c_obj
            !! The pointer to the [[ode_results]] object.
        integer(c_int), intent(in), value :: m
            !! The number of rows in x.
        integer(c_int), intent(in), value :: n
            !! The number of columns in x.
        real(c_double), intent(out) :: x(ldx,n)
            !! The matrix where the results are written.  The matrix is stored
            !! in column-major format with the first column containing the
            !! independent variable values at which the integrator returned
            !! a solution point.  The remaining columns each contain a state
            !! variable.
        integer(c_int), intent(in), value :: ldx
            !! The leading dimension of matrix x.  This value must be at least
            !! m.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        type(ode_results), pointer :: obj

        rst = C_NO_ERROR
        if (.not.c_associated(c_obj)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_pointer(c_obj, obj)
        if (m /= obj%count()) then
            rst = C_MATRIX_ROW_SIZE_ERROR
            return
        end if
        if (n /= obj%state_variable_count() + 1) then
            rst = C_MATRIX_COLUMN_SIZE_ERROR
            return
        end if
        if (ldx < m) then
            rst = C_MATRIX_STRIDE_ERROR
            return
        end if

        x(1:m,1:n) = obj%to_matrix()
    end function

! ------------------------------------------------------------------------------
    function c_write_ode_results(fname, c_obj, delimiter) result(rst) &
        bind(C, name = "c_write_ode_results")
        !! Writes the contents of an [[ode_results]] container to a
        !! delimited text file.
        character(kind = c_char), intent(in) :: fname(*)
            !! The filename.
        type(c_ptr), value, intent(in) :: c_obj
            !! The pointer to the [[ode_results]] object.
        character(kind = c_char), intent(in), value :: delimiter
            !! The delimiter character.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        type(ode_results), pointer :: obj
        character(len = :), allocatable :: filename

        rst = C_NO_ERROR
        if (.not.c_associated(c_obj)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_pointer(c_obj, obj)
        filename = c_to_f_string(fname)
        call write_ode_results(filename, obj, delimiter = delimiter)
    end function

! ------------------------------------------------------------------------------
end module