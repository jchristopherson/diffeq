module diffeq_errors
    !! A collection of routines for handling errors in the DIFFEQ library.
    use iso_fortran_env
    use ferror
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

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
subroutine report_memory_error(err, fcn, flag)
    !! Reports a memory error.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    integer(int32), intent(in) :: flag
        !! The memory status flag.

    ! Local Variables
    character(len = 256) :: msg

    ! Process
    write(msg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error(fcn, trim(msg), DIFFEQ_MEMORY_ALLOCATION_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_matrix_size_error(err, fcn, varname, exprow, expcol, &
    actrow, actcol)
    !! Reports a matrix size error.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    character(len = *), intent(in) :: varname
        !! The offending variable name.
    integer(int32), intent(in) :: exprow
        !! The expected number of rows.
    integer(int32), intent(in) :: expcol
        !! The expected number of columns.
    integer(int32), intent(in) :: actrow
        !! The actual number of rows.
    integer(int32), intent(in) :: actcol
        !! The actual number of columns.

    ! Local Variables
    character(len = 256) :: msg

    ! Process
    write(msg, 100) "The matrix " // varname // " was expected to be (", &
        exprow, "-by-", expcol, "), but was found to be (", actrow, "-by-", &
        actcol, ")."
    call err%report_error(fcn, trim(msg), DIFFEQ_MATRIX_SIZE_ERROR)

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_min_array_size_not_met(err, fcn, varname, minsize, actsize)
    !! Reports an error where the minimum array size was not met.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    character(len = *), intent(in) :: varname
        !! The offending variable name.
    integer(int32), intent(in) :: minsize
        !! The minimum size of the array.
    integer(int32), intent(in) :: actsize
        !! The actual size of the array.

    ! Local Variables
    character(len = 256) :: errmsg

    ! Process
    write(errmsg, 100) "Array " // varname // " must be at least of size ", &
        minsize, ", but was found to be of size ", actsize, "."
    call err%report_error(fcn, trim(errmsg), DIFFEQ_ARRAY_SIZE_ERROR)

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_missing_ode(err, fcn)
    !! Reports a missing ODE routine.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.

    ! Process
    call err%report_error(fcn, "No ODE routine has been supplied.", &
        DIFFEQ_NULL_POINTER_ERROR)
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine report_array_size_error(err, fcn, varname, expsize, actsize)
    !! Reports an array size error.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    character(len = *), intent(in) :: varname
        !! The offending variable name.
    integer(int32), intent(in) :: expsize
        !! The expected size of the array.
    integer(int32), intent(in) :: actsize
        !! The actual size of the array.

    ! Local Variables
    character(len = 256) :: msg

    ! Process
    write(msg, 100) "Array " // varname // " was expected to be of size ", &
        expsize, ", but was found to be of size ", actsize, "."
    call err%report_error(fcn, trim(msg), DIFFEQ_ARRAY_SIZE_ERROR)

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_missing_argument(err, fcn, arg)
    !! Reports a missing argument error.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    character(len = *), intent(in) :: arg
        !! The name of the argument.

    ! Process
    call err%report_error(fcn, "Argument " // arg // " was missing.", &
        DIFFEQ_MISSING_ARGUMENT_ERROR)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_step_size_too_small(err, fcn, x, h)
    !! Reports an error when the step size becomes too small.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which the step size error
        !! occurred.
    real(real64), intent(in) :: h
        !! The step size value.

    ! Process
    character(len = 256) :: msg
    write(msg, 100) "A step size of ", h, " was encountered at x = ", x, &
        ", and is too small to continue integration."
    call err%report_error(fcn, trim(msg), DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR)

    ! Formatting
100 format(A, EN0.3, A, EN0.3, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_excessive_iterations(err, fcn, n, x)
    !! Reports an error when excessive iterations have been made.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    integer(int32), intent(in) :: n
        !! The number of iterations.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which the error occurred.

    ! Process
    character(len = 256) :: msg
    write(msg, 100) "The allowable iteration count was exceeded (iteration ", &
        n, ") at x = ", x, "."
    call err%report_error(fcn, trim(msg), DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)

    ! Formatting
100 format(A, I0, A, EN0.3, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_excessive_integration_steps(err, fcn, n, x)
    !! Reports an error when an excessive  amount integration steps have been 
    !! taken.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    integer(int32), intent(in) :: n
        !! The number of integration steps.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which the error occurred.

    ! Process
    character(len = 256) :: msg
    write(msg, 100) &
        "The allowable number of integration steps was exceeded (step ", &
        n, " at x = ", x, "."
    call err%report_error(fcn, trim(msg), DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)

    ! Formatting
100 format(A, I0, A, EN0.3, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_tolerance_too_small_error(err, fcn)
    !! Reports an error where the requested tolerances are too small.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.

    ! Process
    call err%report_error(fcn, "The requested tolerances are too small.", &
        DIFFEQ_TOLERANCE_TOO_SMALL)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_multiple_convergence_error(err, fcn, n)
    !! Reports an error when multiple convergence tests have failed on a single
    !! integration step.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    integer(int32), intent(in) :: n
        !! The number of failed convergence tests.

    ! Local Variables
    character(len = 256) :: msg

    ! Process
    write(msg, 100) "The integrator has failed to converge ", n, &
        " times on a single integration step."
    call err%report_error(fcn, trim(msg), DIFFEQ_CONVERGENCE_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine report_successive_error_test_failures(err, fcn, n)
    !! Reports an error condition resulting from an integrator failing multiple
    !! successive error tests.
    class(errors), intent(inout) :: err
        !! The error handling object.
    character(len = *), intent(in) :: fcn
        !! The name of the function or subroutine in which the error occurred.
    integer(int32), intent(in) :: n
        !! The number of failed error tests.

    ! Local Variables
    character(len = 256) :: msg

    ! Process
    write(msg, 100) "The integrator has failed ", n, &
        " sequential error tests on a single integration step."
    call err%report_error(fcn, trim(msg), DIFFEQ_CONVERGENCE_ERROR)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end module