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