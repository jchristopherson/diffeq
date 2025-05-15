program test
    use iso_fortran_env
    use diffeq_jacobian_tests
    use diffeq_test_runge_kutta
    use diffeq_test_implicit_rk
    use diffeq_test_vode
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

    rst = test_runge_kutta_45_1()
    if (.not.rst) flag = 5

    rst = test_runge_kutta_45_2()
    if (.not.rst) flag = 6

    rst = test_runge_kutta_45_3()
    if (.not.rst) flag = 7

    rst = test_runge_kutta_23_1()
    if (.not.rst) flag = 8

    rst = test_runge_kutta_23_2()
    if (.not.rst) flag = 9

    rst = test_runge_kutta_23_3()
    if (.not.rst) flag = 10

    rst = test_runge_kutta_853_1()
    if (.not.rst) flag = 11

    rst = test_runge_kutta_853_2()
    if (.not.rst) flag = 12

    rst = test_runge_kutta_853_3()
    if (.not.rst) flag = 13

    rst = test_rosenbrock_1()
    if (.not.rst) flag = 14

    rst = test_rosenbrock_2()
    if (.not.rst) flag = 15

    rst = test_rosenbrock_3()
    if (.not.rst) flag = 16

    rst = test_runge_kutta_with_args()
    if (.not.rst) flag = 17

    rst = test_rosenbrock_with_args()
    if (.not.rst) flag = 18

    rst = test_runge_kutta_dense_with_args()
    if (.not.rst) flag = 19

    rst = test_adams_1()
    if (.not.rst) flag = 20

    ! Output
    stop flag
end program