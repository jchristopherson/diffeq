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

module diffeq_bsrk32_constants
    use iso_fortran_env
    implicit none

    real(real64), parameter :: a21 = 0.5d0
    real(real64), parameter :: a31 = 0.0d0
    real(real64), parameter :: a32 = 0.75d0
    real(real64), parameter :: a41 = 2.0d0 / 9.0d0
    real(real64), parameter :: a42 = 1.0d0 / 3.0d0
    real(real64), parameter :: a43 = 4.0d0 / 9.0d0

    real(real64), parameter :: b1 = 2.0d0 / 9.0d0
    real(real64), parameter :: b2 = 1.0d0 / 3.0d0
    real(real64), parameter :: b3 = 4.0d0 / 9.0d0
    real(real64), parameter :: b4 = 0.0d0

    real(real64), parameter :: b1a = 7.0d0 / 24.0d0
    real(real64), parameter :: b2a = 1.0d0 / 4.0d0
    real(real64), parameter :: b3a = 1.0d0 / 3.0d0
    real(real64), parameter :: b4a = 1.0d0 / 8.0d0

    real(real64), parameter :: c1 = 0.0d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = 0.75d0
    real(real64), parameter :: c4 = 1.0d0

    real(real64), parameter :: e1 = b1a - b1
    real(real64), parameter :: e2 = b2a - b2
    real(real64), parameter :: e3 = b3a - b3
    real(real64), parameter :: e4 = b4a - b4

end module