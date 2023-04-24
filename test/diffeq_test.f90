program test
    use iso_fortran_env
    use diffeq_jacobian_tests
    implicit none

    ! Local Variables
    integer(int32) :: flag
    logical :: rst

    ! Initialization
    flag = 0

    ! Tests
    rst = test_fd_jacobian_1()
    if (.not.rst) flag = 1

    rst = test_fd_jacobian_2()
    if (.not.rst) flag = 2

    ! Output
    stop flag
end program