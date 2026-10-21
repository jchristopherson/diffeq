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

program diffeq_container_io_tests
    use iso_fortran_env, only : real64, iostat_end
    use diffeq, only : ode_results, write_ode_results
    implicit none

    character(len = 32) :: test_name
    type(ode_results) :: results

    call get_command_argument(1, test_name)
    if (len_trim(test_name) == 0) then
        call test_output(",")
        call test_output(";")
        call test_output(achar(9))
        call test_overwrite_and_empty()
    else
        if (trim(test_name) /= "uninitialized") results = ode_results(1)
        select case (trim(test_name))
        case ("quote")
            call write_ode_results("diffeq_writer_quote.txt", results, '"')
        case ("carriage_return")
            call write_ode_results("diffeq_writer_cr.txt", results, achar(13))
        case ("line_feed")
            call write_ode_results("diffeq_writer_lf.txt", results, achar(10))
        case ("uninitialized")
            call write_ode_results("diffeq_writer_uninitialized.txt", results)
        case ("open_error")
            call write_ode_results(".", results)
        case default
            error stop "Unknown writer test"
        end select
        stop 0
    end if

contains
    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len = *), intent(in) :: message
        if (.not.condition) error stop message
    end subroutine

    subroutine test_output(delimiter)
        character, intent(in) :: delimiter
        character(len = *), parameter :: filename = "diffeq_writer_roundtrip.txt"
        type(ode_results) :: container
        real(real64), dimension(3,4) :: expected
        real(real64), dimension(4) :: actual
        real(real64), allocatable, dimension(:,:) :: snapshot
        character(len = 256) :: line
        integer :: unit, flag, index, column, separators

        expected(1,:) = [0.0d0, -1.2345678901234567d0, 1.0d-200, 1.0d200]
        expected(2,:) = [0.125d0, tiny(1.0d0), huge(1.0d0), -tiny(1.0d0)]
        expected(3,:) = [-4.5d0, nearest(1.0d0, 1.0d0), -huge(1.0d0), 0.0d0]
        container = ode_results(3, 7)
        do index = 1, size(expected, 1)
            call container%push(expected(index,1), expected(index,2:))
        end do
        if (delimiter == ",") then
            call write_ode_results(filename, container)
        else
            call write_ode_results(filename, container, delimiter)
        end if

        open(newunit = unit, file = filename, status = "old", action = "read", &
            round = "nearest", iostat = flag)
        call check(flag == 0, "Open writer output")
        do index = 1, size(expected, 1)
            read(unit, "(A)", iostat = flag) line
            call check(flag == 0, "Writer output row count")
            separators = 0
            do column = 1, len_trim(line)
                if (line(column:column) == delimiter) then
                    separators = separators + 1
                    line(column:column) = ","
                end if
            end do
            call check(separators == 3, "Writer delimiter count")
            call check(line(len_trim(line):len_trim(line)) /= ",", "No trailing delimiter")
            call check(scan(line, '"') == 0, "No quotes in numeric output")
            read(line, *, iostat = flag) actual
            call check(flag == 0, "Parse writer numeric output")
            call check(all(actual == expected(index,:)), "Exact double-precision round trip")
        end do
        read(unit, "(A)", iostat = flag) line
        call check(flag == iostat_end, "Writer excludes unused buffer rows")
        close(unit, status = "delete", iostat = flag)
        call check(flag == 0, "Delete writer output")
        call check(container%count() == 3, "Writer preserves count")
        call check(container%get_capacity() == 7, "Writer preserves capacity")
        call check(container%state_variable_count() == 3, "Writer preserves state count")
        snapshot = container%to_matrix()
        call check(all(snapshot == expected), "Writer preserves container contents")
    end subroutine

    subroutine test_overwrite_and_empty()
        character(len = *), parameter :: filename = "diffeq_writer_overwrite.txt"
        type(ode_results) :: container
        character(len = 256) :: line
        real(real64), dimension(2) :: actual
        integer :: unit, flag

        container = ode_results(1, 5)
        call container%push(1.0d0, [2.0d0])
        call container%push(3.0d0, [4.0d0])
        call write_ode_results(filename, container)
        call container%clear_all()
        call container%push(6.0d0, [7.0d0])
        call write_ode_results(filename, container)
        open(newunit = unit, file = filename, status = "old", action = "read", iostat = flag)
        call check(flag == 0, "Open overwritten output")
        read(unit, *, iostat = flag) actual
        call check(flag == 0, "Read single-state output")
        call check(all(actual == [6.0d0, 7.0d0]), "Overwrite replaces contents")
        read(unit, "(A)", iostat = flag) line
        call check(flag == iostat_end, "Overwrite removes old trailing rows")
        close(unit)

        call container%clear_all()
        call write_ode_results(filename, container)
        open(newunit = unit, file = filename, status = "old", action = "read", iostat = flag)
        call check(flag == 0, "Open cleared output")
        read(unit, "(A)", iostat = flag) line
        call check(flag == iostat_end, "Cleared container writes empty file")
        close(unit, status = "delete")

        container = ode_results(2, 3)
        call write_ode_results(filename, container)
        open(newunit = unit, file = filename, status = "old", action = "read", iostat = flag)
        call check(flag == 0, "Open initialized empty output")
        read(unit, "(A)", iostat = flag) line
        call check(flag == iostat_end, "Initialized empty container writes empty file")
        close(unit, status = "delete")
    end subroutine
end program