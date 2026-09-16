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

module diffeq_errors
    !! A collection of routines for handling errors in the DIFFEQ library.
    use iso_fortran_env
    implicit none
    
! ------------------------------------------------------------------------------
    ! Error Flags
    integer(int32), parameter :: DIFFEQ_MEMORY_ALLOCATION_ERROR = 10000
    integer(int32), parameter :: DIFFEQ_NULL_POINTER_ERROR = 10001
    integer(int32), parameter :: DIFFEQ_MATRIX_SIZE_ERROR = 10002
    integer(int32), parameter :: DIFFEQ_ARRAY_SIZE_ERROR = 10003
    integer(int32), parameter :: DIFFEQ_INVALID_INPUT_ERROR = 10004
    integer(int32), parameter :: DIFFEQ_MISSING_ARGUMENT_ERROR = 10005
    integer(int32), parameter :: DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR = 10006
    integer(int32), parameter :: DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR = 10007
    integer(int32), parameter :: DIFFEQ_INVALID_OPERATION_ERROR = 10008
    integer(int32), parameter :: DIFFEQ_TOLERANCE_TOO_SMALL = 10009
    integer(int32), parameter :: DIFFEQ_CONVERGENCE_ERROR = 10010
    integer(int32), parameter :: DIFFEQ_ERROR_TEST_FAILURE = 10011
    integer(int32), parameter :: DIFFEQ_SINGULAR_MATRIX_ERROR = 10012
end module