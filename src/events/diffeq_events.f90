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

module diffeq_events
    use iso_fortran_env
    implicit none
    private
    public :: ode_event
    
    abstract interface
        subroutine ode_event(x, y, f, args)
            !! A callback routine for evaluating event status.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: x
                !! The current value of the independent variable.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element array containing the current state vector.
            real(real64), intent(out), dimension(:) :: f
                !! An N-element array containing the values of the event
                !! functions.  The solver attempts to find the situation 
                !! where f equals zero thereby signaling an event.
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to pass information
                !! in and out of the routine.
        end subroutine
    end interface

contains
! ------------------------------------------------------------------------------
    pure function event_occurred(y1, y2) result(rst)
        !! Tests to see if an event occurred by checking for a change in sign
        !! in any of the state vector components.
        real(real64), intent(in), dimension(:) :: y1
            !! An N-element array containing the state vector at some value x.
        real(real64), intent(in), dimension(:) :: y2
            !! An N-element array containing the state vector at one integration
            !! step away from y1.
        logical :: rst
            !! Returns true if there is a sign change between any of the
            !! components of the two state vectors; else, false, if no sign
            !! change occurred.  A value of 0.0 will also return true.

        rst = any(y1 * y2 <= 0.0d0)
    end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module