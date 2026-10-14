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

module diffeq_containers
    use iso_fortran_env
    use diffeq_errors
    implicit none
    private
    public :: ode_results

    integer, parameter :: BUFFER_CAPACITY = 1000

    type ode_results
        !! A container for ODE results.
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
            !! An M-by-N+1 matrix containing the values of the independent
            !! variable at which the solution is computed in column 1.  The
            !! remaining N columns contain the state variables.  Notice, M
            !! may be greater than the actual number of values stored as 
            !! this is a buffer or values and is typically kept larger than 
            !! actually needed to minimize the number of reallocations.
        integer, private :: m_count = 0
            !! The number of solution points actually stored in the solution
            !! buffers.
    contains
        procedure, public :: state_variable_count => or_get_state_count
        procedure, public :: count => or_get_count
        procedure, public :: set_capacity => or_set_capacity
    end type

    interface ode_results
        module procedure :: or_init_empty
    end interface

contains
! ------------------------------------------------------------------------------
pure function or_init_empty(n) result(rst)
    !! Initializes an empty [[ode_results]] container.
    integer, intent(in) :: n
        !! The number of state variables.  A DIFFEQ_INVALID_INPUT_ERROR is
        !! thrown if this value is less than 1.
    type(ode_results) :: rst
        !! The [[ode_results]] object.
    
    if (n < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    allocate(rst%m_buffer(BUFFER_CAPACITY, n + 1))
    rst%m_count = 0
end function

! ------------------------------------------------------------------------------
pure function or_get_count(this) result(rst)
    !! Gets the number of solution points.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer :: rst
        !! The number of solution points.
    rst = this%m_count
end function

! ------------------------------------------------------------------------------
pure function or_get_state_count(this) result(rst)
    !! Gets the number of state variables.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer :: rst
        !! The number of state variables.
    if (allocated(this%m_buffer)) then
        rst = size(this%m_buffer, 2) - 1
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine or_set_capacity(this, n)
    !! Adjusts the capacity of the container.  If the container is not yet 
    !! initialized a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(inout) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: n
        !! The new capacity.  This value must be greater than one or a
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown.

    ! Local Variables
    integer :: nrows, ncols
    real(real64), allocatable, dimension(:,:) :: temp

    ! Input Checking
    if (n < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR

    ! Process
    call move_alloc(this%m_buffer, temp)
end subroutine

! program demo_move_alloc
!     implicit none
!     integer, allocatable :: a(:), b(:)

!     ! Allocate and initialize A
!     allocate(a(3))
!     a = [1, 2, 3]

!     ! Move allocation from A → B
!     call move_alloc(a, b)

!     ! A is now deallocated; B owns the data
!     print *, "allocated(a) = ", allocated(a)
!     print *, "allocated(b) = ", allocated(b)
!     print *, "b = ", b
! end program demo_move_alloc


! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module