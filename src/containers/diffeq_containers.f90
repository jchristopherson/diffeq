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
    public :: ode_results_statistics
    public :: ode_results
    public :: size

    integer, parameter :: BUFFER_CAPACITY = 1000

    type ode_results_statistics
        integer :: successful_steps = 0
            !! The number of successful steps taken by the integrator.
        integer :: failed_steps = 0
            !! The number of unsuccessful steps taken by the integrator.
        integer :: jacobian_evaluations = 0
            !! The number of Jacobian matrix evaluations.
        integer :: function_evaluations = 0
            !! The number of function evaluations.
    end type

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
        integer, private :: m_state_count = 0
            !! The number of state variables.
        type(ode_results_statistics), public :: statistics
            !! Statistical information regarding the integration process.
    contains
        procedure, public :: state_variable_count => or_get_state_count
        procedure, public :: count => or_get_count
        procedure, public :: get_capacity => or_get_capacity
        procedure, public :: set_capacity => or_set_capacity
        generic, public :: state => or_whole_state, or_state_vector
        procedure, private :: or_whole_state
        procedure, private :: or_state_vector
        procedure, public :: selected_state => or_state_by_index
        procedure, public :: solution_point => or_solution_point_by_index
        procedure, public :: push => or_push_state
        procedure, public :: to_matrix => or_to_matrix
        procedure, public :: clear_all => or_clear_all
        procedure, public :: solution_points => or_solution_points
    end type

    interface ode_results
        module procedure :: or_init_empty
    end interface

    interface size
        module procedure :: or_size
    end interface

contains
! ------------------------------------------------------------------------------
pure function or_size(x, dim) result(rst)
    !! Returns the size of the requested dimension of a [[ode_results]] object.
    class(ode_results), intent(in) :: x
        !! The [[ode_results]] object.
    integer, intent(in) :: dim
        !! The dimension.  A value of 1 returns the number of solution points,
        !! while a value of 2 returns the number of state variables.  Any other
        !! input is disregarded, and a value of 0 is returned.
    integer :: rst
        !! The requested size.

    select case (dim)
    case (1)
        rst = x%count()
    case (2)
        rst = x%state_variable_count()
    case default
        rst = 0
    end select
end function

! ******************************************************************************
! ODE_RESULTS MEMBERS
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
    rst%m_state_count = n
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
        rst = this%m_state_count
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
pure function or_state_by_index(this, index) result(rst)
    !! Gets the requested state variables by index.  If this routine is called
    !! on an uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: index
        !! The one-based index for the state to return.  The index value must
        !! lie within the range [1, count()], else a DIFFEQ_INVALID_INPUT_ERROR
        !! is thrown.
    real(real64), allocatable, dimension(:) :: rst
        !! The state vector.

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (index < 1 .or. index > this%count()) error stop DIFFEQ_INVALID_INPUT_ERROR

    rst = this%m_buffer(index,2:)
end function

! ------------------------------------------------------------------------------
pure function or_solution_point_by_index(this, index) result(rst)
    !! Gets the requested solution point by index.  If this routine is called
    !! on an uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: index
        !! The one-based index for the state to return.  The index value must
        !! lie within the range [1, count()], else a DIFFEQ_INVALID_INPUT_ERROR
        !! is thrown.
    real(real64) :: rst
        !! The solution point.

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (index < 1 .or. index > this%count()) error stop DIFFEQ_INVALID_INPUT_ERROR

    rst = this%m_buffer(index,1)
end function

! ------------------------------------------------------------------------------
pure function or_whole_state(this) result(rst)
    !! Gets the entire contents of the container.  If this routine is called
    !! on an uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An M-by-N element matrix where the state variable values are 
        !! returned.  M is equal to count(), and N is equal to
        !! state_variable_count().

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR    
    rst = this%m_buffer(1:this%count(), 2:this%state_variable_count() + 1)
end function

! ------------------------------------------------------------------------------
pure function or_state_vector(this, index) result(rst)
    !! Gets an entire state vector.   If this routine is called on an 
    !! uninitialized container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    integer, intent(in) :: index
        !! The one-based index for the state vector to return.  The index value 
        !! must lie within the range [1, state_variable_count()], else a 
        !! DIFFEQ_INVALID_INPUT_ERROR is thrown.
    real(real64), allocatable, dimension(:) :: rst
        !! An M-element array where the state vector will be returned.  M is 
        !! equal to count().

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    if (index < 1 .or. index > this%state_variable_count()) error stop DIFFEQ_INVALID_INPUT_ERROR

    rst = this%m_buffer(1:this%m_count,index + 1)
end function

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
        !! number of state variables.  The same applies after clear_all().

    integer :: capacity

    if (.not.allocated(this%m_buffer)) then
        !! This is brand new.  Initialize the collection
        allocate(this%m_buffer(BUFFER_CAPACITY, size(y) + 1))
        this%m_count = 0
        this%m_state_count = size(y)
    else if (this%m_state_count == 0) then
        if (size(this%m_buffer, 2) /= size(y) + 1) then
            capacity = this%get_capacity()
            deallocate(this%m_buffer)
            allocate(this%m_buffer(capacity, size(y) + 1))
        end if
        this%m_state_count = size(y)
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
    rst = this%m_buffer(1:this%count(),1:this%state_variable_count() + 1)
end function

! ------------------------------------------------------------------------------
subroutine or_clear_all(this)
    !! Clears the contents of this container, but does not alter its capacity.
    !! The state variable count is reset; the next push sets the new count.
    class(ode_results), intent(inout) :: this
        !! The [[ode_results]] object.
    this%m_count = 0
    this%m_state_count = 0
end subroutine

! ------------------------------------------------------------------------------
pure function or_solution_points(this) result(rst)
    !! Gets a vector containing the independent variable values at which the
    !! solution was computed.  If this routine is called on an uninitialized 
    !! container a DIFFEQ_INVALID_OPERATION_ERROR is thrown.
    class(ode_results), intent(in) :: this
        !! The [[ode_results]] object.
    real(real64), allocatable, dimension(:) :: rst
        !! An M-element array where the solution points will be written.  M is 
        !! equal to count().

    if (.not.allocated(this%m_buffer)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    rst = this%m_buffer(1:this%count(),1)
end function

! ------------------------------------------------------------------------------
end module