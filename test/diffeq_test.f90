program test
    use iso_fortran_env
    use diffeq_jacobian_tests
    use diffeq_test_runge_kutta
    use diffeq_test_implicit_rk
    use diffeq_test_bdf
    use diffeq_test_pece
    implicit none

    ! Local Variables
    integer(int32) :: flag
    logical :: rst

    ! Initialization
    flag = 0

    ! Tests
    rst = test_state_variable_tolerances()
        if (.not.rst) flag = max(flag, 42)

    rst = test_step_size_limits()
    if (.not.rst) flag = max(flag, 48)

    rst = test_implicit_rk_state_tolerances()
        if (.not.rst) flag = max(flag, 43)

    rst = test_analytical_jacobian_usage()
    if (.not.rst) flag = max(flag, 49)

    rst = test_stiff_vanderpol()
    if (.not.rst) flag = max(flag, 51)

    rst = test_bdf_state_tolerances()
        if (.not.rst) flag = max(flag, 44)

    rst = test_adams_state_tolerances()
        if (.not.rst) flag = max(flag, 45)

    rst = test_reverse_and_overshoot()
        if (.not.rst) flag = max(flag, 46)

    rst = test_bdf_reverse()
        if (.not.rst) flag = max(flag, 47)

    rst = test_tsitouras_54()
        if (.not.rst) flag = max(flag, 1)

    rst = test_tsitouras_54_dense()
        if (.not.rst) flag = max(flag, 2)

    rst = test_fd_jacobian_1()
        if (.not.rst) flag = max(flag, 1)

    rst = test_fd_step_setting()
    if (.not.rst) flag = max(flag, 50)

    rst = test_fd_jacobian_2()
        if (.not.rst) flag = max(flag, 2)

    rst = test_fd_jacobian_3()
        if (.not.rst) flag = max(flag, 3)

    rst = test_fd_jacobian_4()
        if (.not.rst) flag = max(flag, 4)

    rst = test_runge_kutta_45_1()
        if (.not.rst) flag = max(flag, 5)

    rst = test_runge_kutta_45_2()
        if (.not.rst) flag = max(flag, 6)

    rst = test_runge_kutta_45_3()
        if (.not.rst) flag = max(flag, 7)

    rst = test_runge_kutta_23_1()
        if (.not.rst) flag = max(flag, 8)

    rst = test_runge_kutta_23_2()
        if (.not.rst) flag = max(flag, 9)

    rst = test_runge_kutta_23_3()
        if (.not.rst) flag = max(flag, 10)

    rst = test_runge_kutta_853_1()
        if (.not.rst) flag = max(flag, 11)

    rst = test_runge_kutta_853_2()
        if (.not.rst) flag = max(flag, 12)

    rst = test_runge_kutta_853_3()
        if (.not.rst) flag = max(flag, 13)

    rst = test_rosenbrock_1()
        if (.not.rst) flag = max(flag, 14)

    rst = test_rosenbrock_2()
        if (.not.rst) flag = max(flag, 15)

    rst = test_rosenbrock_3()
        if (.not.rst) flag = max(flag, 16)

    rst = test_rosenbrock_mass_matrix()
        if (.not.rst) flag = max(flag, 17)

    rst = test_rosenbrock_with_args()
        if (.not.rst) flag = max(flag, 18)

    rst = test_kennedy_carpenter_4()
        if (.not.rst) flag = max(flag, 19)

    rst = test_kennedy_carpenter_5()
        if (.not.rst) flag = max(flag, 20)

    rst = test_kennedy_carpenter_mass_matrix()
        if (.not.rst) flag = max(flag, 21)

    rst = test_kennedy_carpenter_singular_mass_matrix()
        if (.not.rst) flag = max(flag, 22)

    rst = test_runge_kutta_dense_with_args()
        if (.not.rst) flag = max(flag, 23)

    rst = test_bdf_1()
        if (.not.rst) flag = max(flag, 24)

    rst = test_bdf_2()
        if (.not.rst) flag = max(flag, 25)

    rst = test_bdf_dense()
        if (.not.rst) flag = max(flag, 26)

    rst = test_bdf_all_steps()
        if (.not.rst) flag = max(flag, 27)

    rst = test_bdf_mass_matrix()
        if (.not.rst) flag = max(flag, 28)

    rst = test_bdf_singular_mass_matrix()
        if (.not.rst) flag = max(flag, 29)

    rst = test_bdf_singular_mass_matrix_dense()
        if (.not.rst) flag = max(flag, 30)

    rst = test_bdf_with_args()
        if (.not.rst) flag = max(flag, 31)

    rst = test_bdf_order_range()
        if (.not.rst) flag = max(flag, 32)

    rst = test_adams_1()
        if (.not.rst) flag = max(flag, 33)

    rst = test_adams_2()
        if (.not.rst) flag = max(flag, 34)

    rst = test_adams_dense()
        if (.not.rst) flag = max(flag, 35)

    rst = test_adams_all_steps()
        if (.not.rst) flag = max(flag, 36)

    rst = test_adams_mass_matrix()
        if (.not.rst) flag = max(flag, 37)

    rst = test_adams_with_args()
        if (.not.rst) flag = max(flag, 38)

    rst = test_adams_order_range()
        if (.not.rst) flag = max(flag, 39)

    rst = test_adams_matches_bdf()
        if (.not.rst) flag = max(flag, 40)

    rst = test_adams_high_order()
        if (.not.rst) flag = max(flag, 41)

    ! Output
    stop flag
end program