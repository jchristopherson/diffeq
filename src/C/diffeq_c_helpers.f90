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

module diffeq_c_helpers
    use iso_c_binding
    implicit none
    private
    public :: c_strlen
    public :: c_to_f_string

contains
! ------------------------------------------------------------------------------
    pure function c_strlen(cstr) result(rst)
        !! Returns the length of a C-style null-terminated string.
        character(kind = c_char), intent(in) :: cstr(*)
            !! The C string.
        integer :: rst
            !! The length of the string, not including the null character.

        rst = 0
        do while (cstr(rst + 1) /= C_NULL_CHAR)
            rst = rst + 1
        end do
    end function

! ------------------------------------------------------------------------------
    pure function c_to_f_string(cstr) result(rst)
        !! Converts a null-terminated C-style string to a Fortran string.
        character(kind = c_char), intent(in) :: cstr(*)
            !! The C string.
        character(len = :), allocatable :: rst
            !! The resulting Fortran string.

        integer :: i, n
        n = c_strlen(cstr)
        allocate(character(len = n) :: rst)
        do i = 1, n
            rst(i:i) = cstr(i)
        end do
    end function

! ------------------------------------------------------------------------------
end module