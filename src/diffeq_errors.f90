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
end module