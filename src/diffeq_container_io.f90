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

module diffeq_container_io
    use iso_fortran_env
    use diffeq_containers
    use diffeq_errors, only : DIFFEQ_INVALID_INPUT_ERROR, DIFFEQ_INVALID_OPERATION_ERROR
    implicit none
    private
    public :: write_ode_results

contains
! ------------------------------------------------------------------------------
subroutine write_ode_results(fname, rst, delimiter)
    !! Writes the contents of a [[ode_results]] container to a delimited text
    !! file, replacing any existing file.  Each row contains the independent
    !! variable followed by the state variables, without a header.  An
    !! initialized empty container produces an empty file.  Invalid delimiters
    !! cause DIFFEQ_INVALID_INPUT_ERROR; uninitialized containers and file I/O
    !! failures cause DIFFEQ_INVALID_OPERATION_ERROR.
    character(len = *), intent(in) :: fname
        !! The name of the file to write.
    class(ode_results), intent(in) :: rst
        !! The [[ode_results]] container to write.
    character, intent(in), optional :: delimiter
        !! The field delimiter.  The default is a comma.  The delimiter
        !! must not be a double quote, carriage return, or line feed.

    ! Local Variables
    character :: delim
    character(len = 512) :: message
    integer :: unit, flag, close_flag, index, column
    real(real64), allocatable, dimension(:) :: row

    ! Initialization
    delim = ","
    if (present(delimiter)) delim = delimiter
    if (delim == '"' .or. delim == achar(13) .or. delim == achar(10)) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    if (rst%get_capacity() == 0) error stop DIFFEQ_INVALID_OPERATION_ERROR

    allocate(row(rst%state_variable_count() + 1))
    open(newunit = unit, file = fname, status = "replace", action = "write", &
        form = "formatted", decimal = "point", round = "nearest", &
        iostat = flag, iomsg = message)
    if (flag /= 0) then
        write(error_unit, "(A)") trim(message)
        error stop DIFFEQ_INVALID_OPERATION_ERROR
    end if

    do index = 1, rst%count()
        row(1) = rst%solution_point(index)
        row(2:) = rst%selected_state(index)
        write(unit, "(ES24.16E3,*(A,ES24.16E3))", iostat = flag, iomsg = message) &
            row(1), (delim, row(column), column = 2, size(row))
        if (flag /= 0) then
            close(unit, iostat = close_flag)
            write(error_unit, "(A)") trim(message)
            error stop DIFFEQ_INVALID_OPERATION_ERROR
        end if
    end do

    close(unit, iostat = flag, iomsg = message)
    if (flag /= 0) then
        write(error_unit, "(A)") trim(message)
        error stop DIFFEQ_INVALID_OPERATION_ERROR
    end if
end subroutine

! ------------------------------------------------------------------------------
end module