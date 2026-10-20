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

module diffeq_c_index_reduction
    use iso_c_binding
    use iso_fortran_env
    use diffeq_c_function_pointers
    use diffeq_index_reduction
    use diffeq_c_errors
    use diffeq_base
    implicit none
    private
    public :: c_ggl_container

    type, bind(C) :: c_ggl_container
        type(c_ptr) :: model
        type(c_funptr) :: force
        type(c_funptr) :: mass
        type(c_funptr) :: constraints
        type(c_funptr) :: gradient
        type(c_funptr) :: curvature
    end type

contains
! ------------------------------------------------------------------------------
    function c_create_ggl_container(ncoordinates, nconstraints, force, mass, &
        constraints, gradient, curvature, alpha, beta, c_obj) result(rst) &
        bind(C, name = "c_create_ggl_container")
        !! Constructs a new [[ggl_container]] pointer and returns its address
        !! to the supplied C pointer.  Be sure to call c_free_ggl_container
        !! once finished with this object.
        integer(c_int), intent(in), value :: ncoordinates
            !! The number of coordinates (state variables of the ODE component).
            !! This value must be positive-valued; else, a C_INVALID_INPUT_ERROR
            !! is returned.
        integer(c_int), intent(in), value :: nconstraints
            !! The number of constraint equations.  This value must be 
            !! positive-valued; else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), intent(in), value :: force
            !! A pointer to the right-hand-side force routine.  A 
            !! C_NULL_POINTER_ERROR error is returned if NULL.
        type(c_funptr), intent(in), value :: mass
            !! A pointer to the physical mass routine.  Supply NULL in the event
            !! that an identity matrix is desired.
        type(c_funptr), intent(in), value :: constraints
            !! A pointer to the routine containing the algebraic constraints.  A 
            !! C_NULL_POINTER_ERROR error is returned if NULL.
        type(c_funptr), intent(in), value :: gradient
            !! A pointer to the routine that computes the constraint gradient.
            !! A C_NULL_POINTER_ERROR error is returned if NULL.
        type(c_funptr), intent(in), value :: curvature
            !! A pointer to the routine that computes the curvature of the
            !! constraint equations.  A C_NULL_POINTER_ERROR error is returned 
            !! if NULL.
        real(c_double), intent(in), value :: alpha
            !! A non-negative, finite position decay rate.  Zero disables the
            !! damping term without disabling the differentiated constraint.
        real(c_double), intent(in), value :: beta
            !! A non-negative, finite velocity decay rate.  Zero disables the
            !! damping term without disabling the differentiated constraint.
        type(c_ggl_container), intent(out) :: c_obj
            !! A pointer to the newly initialized [[c_ggl_container]].
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.  

        type(ggl_container), pointer :: obj

        rst = C_NO_ERROR
        if (nconstraints < 1) then
            rst = C_INVALID_INPUT_ERROR
            return
        end if
        if (ncoordinates < 1) then
            rst = C_INVALID_INPUT_ERROR
            return
        end if

        allocate(obj)
        call obj%initialize(ncoordinates, nconstraints, alpha = alpha, beta = beta)
        c_obj%model = c_loc(obj)
        
        if (.not.c_associated(force)) then
            rst = C_NULL_POINTER_ERROR
            deallocate(obj)
            return
        end if
        c_obj%force = force

        ! MASS Can be NULL
        c_obj%mass = mass

        if (.not.c_associated(constraints)) then
            rst = C_NULL_POINTER_ERROR
            deallocate(obj)
            return
        end if
        c_obj%constraints = constraints

        if (.not.c_associated(gradient)) then
            rst = C_NULL_POINTER_ERROR
            deallocate(obj)
            return
        end if
        c_obj%gradient = gradient

        if (.not.c_associated(curvature)) then
            rst = C_NULL_POINTER_ERROR
            deallocate(obj)
            return
        end if
        c_obj%curvature = curvature
    end function

! ------------------------------------------------------------------------------
    subroutine c_free_ggl_results(c_obj) bind(C, name = "c_free_ggl_results")
        !! Cleans up resources held by a [[c_ggl_container]] object.
        type(c_ggl_container), intent(inout) :: c_obj
            !! A pointer to the [[c_ggl_container]] object to free.

        type(ggl_container), pointer :: obj
        if (.not.c_associated(c_obj%model)) return
        call c_f_pointer(c_obj%model, obj)
        if (associated(obj)) deallocate(obj)
        c_obj%model = c_null_ptr
        c_obj%constraints = c_null_funptr
        c_obj%curvature = c_null_funptr
        c_obj%force = c_null_funptr
        c_obj%gradient = c_null_funptr
        c_obj%mass = c_null_funptr
    end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module