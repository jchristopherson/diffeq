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

module c_plot_helper
    use iso_c_binding
    use iso_fortran_env
    use fplot_core
    use diffeq_c_helpers
    implicit none

contains
subroutine c_plot_xy(n, x, y, xlbl, ylbl) bind(C, name = "c_plot_xy")
    integer(c_int), intent(in), value :: n
    real(real64), intent(in) :: x(n)
    real(real64), intent(in) :: y(n)
    character(kind = c_char), intent(in) :: xlbl(*)
    character(kind = c_char), intent(in) :: ylbl(*)

    type(plot_2d) :: plt

    call plt%initialize()
    call plt%set_x_axis_title(c_to_f_string(xlbl))
    call plt%set_y_axis_title(c_to_f_string(ylbl))
    call plt%push(x, y)
    call plt%draw()
end subroutine

subroutine c_plot_xy2(n1, n2, x1, y1, name1, x2, y2, name2, xlbl, ylbl) &
    bind(C, name = "c_plot_xy2")
    integer(c_int), intent(in), value :: n1
    integer(c_int), intent(in), value :: n2
    real(c_double), intent(in) :: x1(n1)
    real(c_double), intent(in) :: y1(n1)
    character(kind = c_char), intent(in) :: name1(*)
    real(c_double), intent(in) :: x2(n2)
    real(c_double), intent(in) :: y2(n2)
    character(kind = c_char), intent(in) :: name2(*)
    character(kind = c_char), intent(in) :: xlbl(*)
    character(kind = c_char), intent(in) :: ylbl(*)

    type(plot_2d) :: plt
    class(legend), pointer :: lgnd

    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)

    call plt%set_x_axis_title(c_to_f_string(xlbl))
    call plt%set_y_axis_title(c_to_f_string(ylbl))

    call plt%push(x1, y1, name = c_to_f_string(name1))
    call plt%push(x2, y2, name = c_to_f_string(name2))

    call plt%draw()
end subroutine

end module