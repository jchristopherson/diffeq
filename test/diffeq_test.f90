program test
    use iso_fortran_env
    use diffeq_jacobian_tests
    use diffeq_rkfixed_tests
    use diffeq_expfixed_tests
    use diffeq_adamsfixed_tests
    use diffeq_dprk45_tests
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

    rst = test_fd_jacobian_3()
    if (.not.rst) flag = 3

    rst = test_fd_jacobian_4()
    if (.not.rst) flag = 4

    rst = test_rk4_1()
    if (.not.rst) flag = 5

    rst = test_rk4_2()
    if (.not.rst) flag = 6

    rst = test_expf_1()
    if (.not.rst) flag = 7

    rst = test_expf_2()
    if (.not.rst) flag = 8

    rst = test_adams_1()
    if (.not.rst) flag = 9

    rst = test_adams_2()
    if (.not.rst) flag = 10

    rst = test_dprk45_attempt_step()
    if (.not.rst) flag = 11

    rst = test_dprk45_step()
    if (.not.rst) flag = 12

    rst = test_dprk45_1()
    if (.not.rst) flag = 13

    rst = test_dprk45_2()
    if (.not.rst) flag = 14

    ! Output
    stop flag
end program