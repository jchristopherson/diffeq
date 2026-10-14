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
        procedure, public :: get_capacity => or_get_capacity
        procedure, public :: set_capacity => or_set_capacity
        generic, public :: state => or_state_by_index, or_whole_state, &
            or_state_vector
        procedure, private :: or_state_by_index
        procedure, private :: or_whole_state
        procedure, private :: or_state_vector
        procedure, public :: push => or_push_state
        procedure, public :: to_matrix => or_to_matrix
    end type

    interface ode_results
        module procedure :: or_init_empty
    end interface

contains
! ------------------------------------------------------------------------------
pure function or_init_empty(n, capacity) result(rst)
    !! Initializes an empty [[ode_results]] container.
    integer, intent(in) :: n
        !! The number of state variables.  A DIFFEQ_INVALID_INPUT_ERROR is
        !! thrown if this value is less than 1.
    integer, intent(in), optional :: capacity
        !! An optional input that sets the initial buffer capacity.  If 
        !! supplied, this value must be at least one or a 
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown.
    type(ode_results) :: rst
        !! The [[ode_results]] object.
    
    if (n < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (present(capacity)) then
        if (capacity < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
        allocate(rst%m_buffer(capacity, n + 1))
    else
        allocate(rst%m_buffer(BUFFER_CAPACITY, n + 1))
    end if
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
pure function or_get_capacity(this) result(rst)
    !! Gets the current capacity of the container.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer :: rst
        !! The current capacity.
    if (allocated(this%m_buffer)) then
        rst = size(this%m_buffer, 1)
    else
        rst = 0
    end if
end function

! ------------------------------------------------------------------------------
subroutine or_set_capacity(this, n)
    !! Adjusts the capacity of the container.  If the container is not yet 
    !! initialized a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    !! When shrinking below count(), only the first n solution points are kept.
    class(ode_results), intent(inout) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: n
        !! The new capacity.  This value must be at least one or a
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown.

    ! Local Variables
    integer :: nrows, ncols
    real(real64), allocatable, dimension(:,:) :: temp

    ! Input Checking
    if (n < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR

    ! Process
    if (n == this%get_capacity()) return
    nrows = min(n, this%m_count)
    ncols = size(this%m_buffer, 2)
    allocate(temp(n, ncols))
    temp(1:nrows,:) = this%m_buffer(1:nrows,:)
    call move_alloc(temp, this%m_buffer)
    this%m_count = nrows
end subroutine


! ------------------------------------------------------------------------------
pure subroutine or_state_by_index(this, index, x, y)
    !! Gets the requested state variables by index.  If this routine is called
    !! on an uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: index
        !! The one-based index for the state to return.  The index value must
        !! lie within the range [1, count()], else a DIFFEQ_INVALID_INPUT_ERROR
        !! is thrown.
    real(real64), intent(out) :: x
        !! The independent variable value.
    real(real64), intent(out), dimension(:) :: y
        !! The state vector.  This vector must be correctly sized to accept the
        !! state vector, else a DIFFEQ_ARRAY_SIZE_ERROR is thrown.

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (size(y) /= this%state_variable_count()) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (index < 1 .or. index > this%count()) error stop DIFFEQ_INVALID_INPUT_ERROR

    x = this%m_buffer(index,1)
    y = this%m_buffer(index,2:)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine or_whole_state(this, x, y)
    !! Gets the entire contents of the container.  If this routine is called
    !! on an uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    real(real64), intent(out), dimension(:) :: x
        !! An M-element array where the independent variable values at which the
        !! solution was computed are returned.  M must be equal to count() or
        !! a DIFFEQ_ARRAY_SIZE_ERROR is thrown.
    real(real64), intent(out), dimension(:,:) :: y
        !! An M-by-N element matrix where the state variable values are 
        !! returned.  M must be equal to count(), and N must be equal to
        !! state_variable_count() or a DIFFEQ_MATRIX_SIZE_ERROR is thrown.

    integer :: i

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (size(x) /= this%count()) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (size(y, 1) /= this%count() .or. &
        size(y, 2) /= this%state_variable_count()) &
        error stop DIFFEQ_MATRIX_SIZE_ERROR

    do i = 1, this%count()
        x(i) = this%m_buffer(i,1)
        y(i,:) = this%m_buffer(i,2:)
    end do
end subroutine

! ------------------------------------------------------------------------------
pure subroutine or_state_vector(this, index, yi)
    !! Gets an entire state vector.   If this routine is called on an 
    !! uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: index
        !! The one-based index for the state vector to return.  The index value 
        !! must lie within the range [1, state_variable_count()], else a 
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown.
    real(real64), intent(out), dimension(:) :: yi
        !! An N-element array where the state vector will be returned.  If N
        !! is not equal to count(), a DIFFEQ_ARRAY_SIZE_ERROR will be thrown.

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (index < 1 .or. index > this%state_variable_count()) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (size(yi) /= this%count()) error stop DIFFEQ_ARRAY_SIZE_ERROR

    yi = this%m_buffer(1:this%m_count,index + 1)
end subroutine

! ------------------------------------------------------------------------------
subroutine or_push_state(this, x, y)
    !! Pushes a new state vector into the back of the container.
    class(ode_results), intent(inout) :: this
        !! The [[ode_results]] object.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! The current state vector.  The dimensions of this vector must match
        !! the value returned by state_variable_count(); else, a 
        !! DIFFEQ_ARRAY_SIZE_ERROR will be thrown.  If the container is not
        !! yet initialized, the size of this vector will set the expected 
        !! number of state variables.

    if (.not.allocated(this%m_buffer)) then
        !! This is brand new.  Initialize the collection
        allocate(this%m_buffer(BUFFER_CAPACITY, size(y) + 1))
        this%m_count = 0
    else
        if (size(y) /= this%state_variable_count()) error stop DIFFEQ_ARRAY_SIZE_ERROR
    end if
    if (this%m_count == this%get_capacity()) then
        call this%set_capacity(this%get_capacity() + BUFFER_CAPACITY)
    end if
    this%m_count = this%m_count + 1
    this%m_buffer(this%m_count, 1) = x
    this%m_buffer(this%m_count, 2:) = y
end subroutine

! ------------------------------------------------------------------------------
pure function or_to_matrix(this) result(rst)
    !! Returns a copy of the contents of this container as a matrix with the 
    !! first column containing the independent variable values at which the
    !! solution was computed and the remaining columns containing the state
    !! vectors.  If this routine is called on an uninitialized container a 
    !! DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! The matrix.

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    rst = this%m_buffer(1:this%count(),:)
end function

! ------------------------------------------------------------------------------
end module