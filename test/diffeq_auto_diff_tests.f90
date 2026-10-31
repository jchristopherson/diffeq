program diffeq_auto_diff_tests
    !! Regression tests for values, directional derivatives, and seeded Jacobians.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use diffeq_auto_diff
    implicit none

    call test_assignment()
    call test_arithmetic()
    call test_comparisons()
    call test_intrinsics()
    call test_arrays()
    call test_jacobian()
    print *, "Automatic differentiation tests passed."

contains

    subroutine check(condition, label)
        !! Stop with a descriptive label if a logical expectation fails.
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not.condition) then
            print *, "FAIL: ", label
            error stop 1
        end if
    end subroutine

    subroutine check_real(actual, expected, label)
        !! Check a finite real result against a scaled analytic reference.
        real(real64), intent(in) :: actual, expected
        character(*), intent(in) :: label
        call check(ieee_is_finite(actual), label // " finite")
        if (abs(actual - expected) > 5.0d-13 * max(1.0d0, abs(expected))) then
            print *, "FAIL: ", label, " actual: ", actual, " expected: ", expected
            error stop 1
        end if
    end subroutine

    subroutine check_ad(actual, value, derivative, label)
        !! Check both components, including nonunit and negative derivative seeds.
        type(real_ad), intent(in) :: actual
        real(real64), intent(in) :: value, derivative
        character(*), intent(in) :: label
        call check_real(actual%x, value, label // " value")
        call check_real(actual%dx, derivative, label // " derivative")
    end subroutine

    subroutine test_assignment()
        !! Verify copying, constant seed clearing, extraction, and scalar expansion.
        type(real_ad) :: a, b(3)
        real(real64) :: values(3), scalar

        a = real_ad(2.0d0, -3.0d0)
        b = a
        call check(all(b%x == 2.0d0) .and. all(b%dx == -3.0d0), "AD broadcast assignment")
        a = 4.0d0
        call check_ad(a, 4.0d0, 0.0d0, "real assignment clears seed")
        a%dx = 7.0d0
        a = 2
        call check_ad(a, 2.0d0, 0.0d0, "integer assignment clears seed")
        b = [1.0d0, 2.0d0, 3.0d0]
        call check(all(b%dx == 0.0d0), "real array assignment clears seeds")
        b%dx = -5.0d0
        values = b
        call check(all(values == [1.0d0, 2.0d0, 3.0d0]), "real array extraction")
        b = [3, 2, 1]
        call check(all(b%x == [3.0d0, 2.0d0, 1.0d0]) .and. all(b%dx == 0.0d0), &
            "integer array assignment")
        a = real_ad(3.0d0, 9.0d0)
        scalar = a
        call check_real(scalar, 3.0d0, "scalar extraction")
    end subroutine

    subroutine test_arithmetic()
        !! Exercise every arithmetic signature and sensitive power-rule cases.
        type(real_ad) :: a, b, z, negative

        a = real_ad(2.0d0, 3.0d0)
        b = real_ad(4.0d0, -2.0d0)
        call check_ad(+a, 2.0d0, 3.0d0, "unary plus")
        call check_ad(-a, -2.0d0, -3.0d0, "unary minus")
        call check_ad(a + b, 6.0d0, 1.0d0, "AD addition")
        call check_ad(a - b, -2.0d0, 5.0d0, "AD subtraction")
        call check_ad(a * b, 8.0d0, 8.0d0, "AD multiplication")
        call check_ad(a / b, 0.5d0, 1.0d0, "AD division")
        call check_ad(a + 4.0d0, 6.0d0, 3.0d0, "AD + real")
        call check_ad(4.0d0 + a, 6.0d0, 3.0d0, "real + AD")
        call check_ad(a - 4.0d0, -2.0d0, 3.0d0, "AD - real")
        call check_ad(4.0d0 - a, 2.0d0, -3.0d0, "real - AD")
        call check_ad(a * 4.0d0, 8.0d0, 12.0d0, "AD * real")
        call check_ad(4.0d0 * a, 8.0d0, 12.0d0, "real * AD")
        call check_ad(a / 4.0d0, 0.5d0, 0.75d0, "AD / real")
        call check_ad(4.0d0 / a, 2.0d0, -3.0d0, "real / AD")
        call check_ad(a + 4, 6.0d0, 3.0d0, "AD + integer")
        call check_ad(4 + a, 6.0d0, 3.0d0, "integer + AD")
        call check_ad(a - 4, -2.0d0, 3.0d0, "AD - integer")
        call check_ad(4 - a, 2.0d0, -3.0d0, "integer - AD")
        call check_ad(a * 4, 8.0d0, 12.0d0, "AD * integer")
        call check_ad(4 * a, 8.0d0, 12.0d0, "integer * AD")
        call check_ad(a / 4, 0.5d0, 0.75d0, "AD / integer")
        call check_ad(4 / a, 2.0d0, -3.0d0, "integer / AD")
        call check_ad(a ** b, 16.0d0, 96.0d0 - 32.0d0 * log(2.0d0), "AD ** AD")
        call check_ad(a ** 2.5d0, 2.0d0**2.5d0, 7.5d0 * 2.0d0**1.5d0, "AD ** real")
        call check_ad(4.0d0 ** a, 16.0d0, 48.0d0 * log(4.0d0), "real ** AD")
        call check_ad(a ** 3, 8.0d0, 36.0d0, "AD ** integer")
        call check_ad(4 ** a, 16.0d0, 48.0d0 * log(4.0d0), "integer ** AD")
        call check_ad(a ** (-2), 0.25d0, -0.75d0, "negative integer power")
        call check_ad(a ** (-2.0d0), 0.25d0, -0.75d0, "negative real power")

        z = real_ad(0.0d0, 3.0d0)
        call check_ad(z ** 0, 1.0d0, 0.0d0, "zero ** integer zero")
        call check_ad(z ** 0.0d0, 1.0d0, 0.0d0, "zero ** real zero")
        call check_ad(z ** 1, 0.0d0, 3.0d0, "zero ** integer one")
        call check_ad(z ** 1.0d0, 0.0d0, 3.0d0, "zero ** real one")
        call check_ad(z ** 2, 0.0d0, 0.0d0, "zero ** integer two")
        call check_ad(z ** 2.0d0, 0.0d0, 0.0d0, "zero ** real two")
        call check_ad(z ** a, 0.0d0, 0.0d0, "zero AD base with variable exponent")
        call check_ad(0.0d0 ** a, 0.0d0, 0.0d0, "zero real base with variable exponent")
        call check_ad(0 ** a, 0.0d0, 0.0d0, "zero integer base with variable exponent")
        z%dx = 0.0d0
        call check_ad(z ** 0.5d0, 0.0d0, 0.0d0, "inactive square-root power")
        call check_ad(a ** z, 1.0d0, 0.0d0, "AD exponent zero")

        negative = real_ad(-2.0d0, 3.0d0)
        call check_ad(negative ** 3, -8.0d0, 36.0d0, "negative base integer power")
        call check_ad(negative ** (-2), 0.25d0, 0.75d0, "negative base reciprocal power")
        call check_ad(negative ** 3.0d0, -8.0d0, 36.0d0, "negative base real integral power")
        b = 3.0d0
        call check_ad(negative ** b, -8.0d0, 36.0d0, "negative base constant AD exponent")
        call check_ad((-2.0d0) ** b, -8.0d0, 0.0d0, "negative real base inactive exponent")
        call check_ad(a / a, 1.0d0, 0.0d0, "quotient cancellation")
        call check_ad(a - a, 0.0d0, 0.0d0, "difference cancellation")
    end subroutine

    subroutine test_comparisons()
        !! Compare all operand pairings below, at, and above equality.
        type(real_ad) :: a, b
        logical :: expected(6)
        integer :: i

        a = real_ad(2.0d0, 3.0d0)
        do i = 1, 3
            b = real_ad(real(i, real64), -7.0d0)
            expected = [2.0d0 == b%x, 2.0d0 /= b%x, 2.0d0 < b%x, &
                2.0d0 <= b%x, 2.0d0 > b%x, 2.0d0 >= b%x]
            call check(all([a == b, a /= b, a < b, a <= b, a > b, a >= b] .eqv. expected), &
                "AD comparisons ignore derivatives")
            call check(all([a == b%x, a /= b%x, a < b%x, a <= b%x, a > b%x, a >= b%x] &
                .eqv. expected), "AD-real comparisons")
            call check(all([2.0d0 == b, 2.0d0 /= b, 2.0d0 < b, 2.0d0 <= b, 2.0d0 > b, 2.0d0 >= b] &
                .eqv. expected), "real-AD comparisons")
            call check(all([a == i, a /= i, a < i, a <= i, a > i, a >= i] .eqv. expected), &
                "AD-integer comparisons")
            call check(all([2 == b, 2 /= b, 2 < b, 2 <= b, 2 > b, 2 >= b] .eqv. expected), &
                "integer-AD comparisons")
        end do
        call check(a .eq. 2 .and. a .ne. 3 .and. a .lt. 3 .and. a .le. 2 .and. &
            a .gt. 1 .and. a .ge. 2, "dotted comparison aliases")
    end subroutine

    subroutine test_intrinsics()
        !! Check each intrinsic against its analytic derivative and branch convention.
        type(real_ad) :: a, b, c
        real(real64) :: x, d, radius2, factor
        integer :: sx, sy

        x = 0.4d0
        d = -2.5d0
        a = real_ad(x, d)
        call check_ad(abs(a), x, d, "abs positive")
        call check_ad(abs(-a), x, d, "abs negative")
        call check_ad(abs(real_ad(0.0d0, d)), 0.0d0, 0.0d0, "abs zero convention")
        call check_ad(sqrt(a), sqrt(x), d / (2.0d0 * sqrt(x)), "sqrt")
        call check_ad(exp(a), exp(x), d * exp(x), "exp")
        call check_ad(log(a), log(x), d / x, "log")
        call check_ad(log10(a), log10(x), d / (x * log(10.0d0)), "log10")
        call check_ad(sin(a), sin(x), d * cos(x), "sin")
        call check_ad(cos(a), cos(x), -d * sin(x), "cos")
        call check_ad(tan(a), tan(x), d / cos(x)**2, "tan")
        call check_ad(asin(a), asin(x), d / sqrt(1.0d0 - x*x), "asin")
        call check_ad(acos(a), acos(x), -d / sqrt(1.0d0 - x*x), "acos")
        call check_ad(atan(a), atan(x), d / (1.0d0 + x*x), "atan")
        call check_ad(sinh(a), sinh(x), d * cosh(x), "sinh")
        call check_ad(cosh(a), cosh(x), d * sinh(x), "cosh")
        call check_ad(tanh(a), tanh(x), d / cosh(x)**2, "tanh")
        call check_ad(tanh(real_ad(1000.0d0, d)), 1.0d0, 0.0d0, "large positive tanh")
        call check_ad(tanh(real_ad(-1000.0d0, d)), -1.0d0, 0.0d0, "large negative tanh")
        c = tanh(real_ad(20.0d0, 1.0d0))
        call check_real(c%dx * cosh(20.0d0)**2, 1.0d0, "saturated tanh retains derivative")
        call check_ad(asinh(a), asinh(x), d / sqrt(1.0d0 + x*x), "asinh")
        call check_ad(atanh(a), atanh(x), d / (1.0d0 - x*x), "atanh")
        factor = 2.0d0 / sqrt(acos(-1.0d0)) * exp(-x*x)
        call check_ad(erf(a), erf(x), d * factor, "erf")
        call check_ad(erfc(a), erfc(x), -d * factor, "erfc")
        call check_ad(acosh(a + 1), acosh(x + 1.0d0), d / sqrt((x + 1.0d0)**2 - 1.0d0), "acosh")
        call check_ad(sqrt(real_ad(0.0d0, 0.0d0)), 0.0d0, 0.0d0, "inactive sqrt at zero")
        call check_ad(asin(real_ad(1.0d0, 0.0d0)), asin(1.0d0), 0.0d0, "inactive asin endpoint")
        call check_ad(acos(real_ad(-1.0d0, 0.0d0)), acos(-1.0d0), 0.0d0, "inactive acos endpoint")
        call check_ad(acosh(real_ad(1.0d0, 0.0d0)), 0.0d0, 0.0d0, "inactive acosh endpoint")

        do sx = -1, 1, 2
            do sy = -1, 1, 2
                a = real_ad(sy * 2.0d0, -3.0d0)
                b = real_ad(sx * 3.0d0, 4.0d0)
                radius2 = a%x**2 + b%x**2
                call check_ad(atan2(a, b), atan2(a%x, b%x), &
                    (b%x * a%dx - a%x * b%dx) / radius2, "atan2 quadrants")
                call check_ad(atan2(a, b%x), atan2(a%x, b%x), b%x * a%dx / radius2, "atan2 AD-real")
                call check_ad(atan2(a%x, b), atan2(a%x, b%x), -a%x * b%dx / radius2, "atan2 real-AD")
                c = atan2(a, b)
                call check_ad(atan(a, b), c%x, c%dx, "two-argument atan AD-AD")
                c = atan2(a, b%x)
                call check_ad(atan(a, b%x), c%x, c%dx, "two-argument atan AD-real")
                c = atan2(a%x, b)
                call check_ad(atan(a%x, b), c%x, c%dx, "two-argument atan real-AD")
            end do
        end do
        call check_ad(atan2(real_ad(0.0d0, 2.0d0), 3.0d0), 0.0d0, 2.0d0 / 3.0d0, "atan2 y zero")
        call check_ad(atan2(3.0d0, real_ad(0.0d0, 2.0d0)), atan2(3.0d0, 0.0d0), &
            -2.0d0 / 3.0d0, "atan2 x zero")

        a = real_ad(2.0d0, 3.0d0)
        b = real_ad(4.0d0, -5.0d0)
        call check_ad(min(a, b), 2.0d0, 3.0d0, "min first")
        call check_ad(min(b, a), 2.0d0, 3.0d0, "min second")
        call check_ad(max(a, b), 4.0d0, -5.0d0, "max second")
        call check_ad(max(b, a), 4.0d0, -5.0d0, "max first")
        call check_ad(min(a, 1.0d0), 1.0d0, 0.0d0, "min AD-real constant")
        call check_ad(min(a, 3.0d0), 2.0d0, 3.0d0, "min AD-real active")
        call check_ad(min(1.0d0, a), 1.0d0, 0.0d0, "min real-AD constant")
        call check_ad(min(3.0d0, a), 2.0d0, 3.0d0, "min real-AD active")
        call check_ad(max(a, 3.0d0), 3.0d0, 0.0d0, "max AD-real constant")
        call check_ad(max(a, 1.0d0), 2.0d0, 3.0d0, "max AD-real active")
        call check_ad(max(3.0d0, a), 3.0d0, 0.0d0, "max real-AD constant")
        call check_ad(max(1.0d0, a), 2.0d0, 3.0d0, "max real-AD active")
        b%x = a%x
        call check_ad(min(a, b), 2.0d0, 3.0d0, "min tie first seed")
        call check_ad(max(a, b), 2.0d0, 3.0d0, "max tie first seed")
        call check_ad(min(a, 2.0d0), 2.0d0, 3.0d0, "min AD-real tie")
        call check_ad(max(a, 2.0d0), 2.0d0, 3.0d0, "max AD-real tie")
        call check_ad(min(2.0d0, a), 2.0d0, 0.0d0, "min real-AD tie")
        call check_ad(max(2.0d0, a), 2.0d0, 0.0d0, "max real-AD tie")

        call check_ad(sign(a, b), 2.0d0, 3.0d0, "sign positive")
        call check_ad(sign(a, -b), -2.0d0, -3.0d0, "sign negative")
        call check_ad(sign(-a, -b), -2.0d0, -3.0d0, "sign negative magnitude")
        call check_ad(sign(a, -1.0d0), -2.0d0, -3.0d0, "sign AD-real")
        call check_ad(sign(-2.0d0, b), 2.0d0, 0.0d0, "sign real-AD ignores seed")
        call check_ad(sign(a, real_ad(0.0d0, 8.0d0)), 2.0d0, 3.0d0, "sign at positive zero")
        call check_ad(sign(real_ad(0.0d0, 8.0d0), -b), -0.0d0, 0.0d0, "sign zero magnitude")

        ! Large arguments should not overflow intermediate squared denominators.
        a = real_ad(1.0d200, 1.0d200)
        call check_ad(asinh(a), asinh(a%x), 1.0d0, "large asinh")
        call check_ad(acosh(a), acosh(a%x), 1.0d0, "large acosh")
        c = atan(a)
        call check(c%dx > 0.0d0, "large atan derivative not lost")
        call check_real(c%dx / 1.0d-200, 1.0d0, "large atan derivative")
        b = real_ad(1.0d200, -1.0d200)
        call check_ad(atan2(a, b), atan2(a%x, b%x), 1.0d0, "large atan2")
        a = real_ad(1.0d308, 1.0d308)
        c = atan(a)
        call check(c%dx > 0.0d0, "very large atan derivative not lost")
        call check_real(c%dx / 1.0d-308, 1.0d0, "very large atan derivative")
    end subroutine

    subroutine test_arrays()
        !! Verify elemental intrinsic calls, arithmetic, comparisons, and rank two.
        type(real_ad) :: a(3), b(3), matrix(2,2), result(2,2)
        real(real64) :: x(3)
        integer :: i

        x = [0.2d0, 0.4d0, 0.6d0]
        a = x
        a%dx = [1.0d0, -2.0d0, 3.0d0]
        b = sin(a) + a**2 / 2
        do i = 1, 3
            call check_ad(b(i), sin(x(i)) + x(i)**2 / 2.0d0, &
                (cos(x(i)) + x(i)) * a(i)%dx, "elemental expression")
        end do
        call check(all(a < 1.0d0) .and. all(0 < a) .and. all(a == x), "array comparisons")
        b = max(a, 0.3d0)
        call check_ad(b(1), 0.3d0, 0.0d0, "elemental branch constant")
        call check_ad(b(2), 0.4d0, -2.0d0, "elemental branch active")
        matrix = 2.0d0
        matrix%dx = 3.0d0
        result = exp(matrix) * matrix
        call check(all(abs(result%x - 2.0d0 * exp(2.0d0)) < 1.0d-12), "rank-two values")
        call check(all(abs(result%dx - 9.0d0 * exp(2.0d0)) < 1.0d-12), "rank-two derivatives")
    end subroutine

    pure function rhs(y) result(f)
        !! Compute a nonlinear RHS using only AD operations and a value-based branch.
        type(real_ad), intent(in) :: y(2)
        type(real_ad) :: f(2)
        f(1) = y(1)**2 * y(2) + sin(y(2))
        if (y(1) > 0) then
            f(2) = exp(y(1)) / y(2) + y(1)**y(2)
        else
            f(2) = -y(1)
        end if
    end function

    subroutine test_jacobian()
        !! Recover both analytic Jacobian columns and a general directional product.
        type(real_ad) :: y(2), f(2)
        real(real64) :: point(2), jac(2,2), exact(2,2), direction(2)
        integer :: col, row

        point = [1.2d0, 0.7d0]
        exact(1,1) = 2.0d0 * point(1) * point(2)
        exact(1,2) = point(1)**2 + cos(point(2))
        exact(2,1) = exp(point(1)) / point(2) + point(2) * point(1)**(point(2) - 1.0d0)
        exact(2,2) = -exp(point(1)) / point(2)**2 + point(1)**point(2) * log(point(1))
        do col = 1, 2
            y = point
            y(col)%dx = 1.0d0
            f = rhs(y)
            jac(:,col) = f%dx
            call check_real(f(1)%x, point(1)**2 * point(2) + sin(point(2)), "Jacobian RHS value 1")
            call check_real(f(2)%x, exp(point(1)) / point(2) + point(1)**point(2), "Jacobian RHS value 2")
        end do
        do col = 1, 2
            do row = 1, 2
                call check_real(jac(row,col), exact(row,col), "analytic Jacobian entry")
            end do
        end do
        direction = [2.0d0, -3.0d0]
        y = point
        y%dx = direction
        f = rhs(y)
        do row = 1, 2
            call check_real(f(row)%dx, dot_product(exact(row,:), direction), "Jacobian-vector product")
        end do
        y(1) = real_ad(-1.0d0, 2.0d0)
        f = rhs(y)
        call check_ad(f(2), 1.0d0, -2.0d0, "negative RHS branch")
    end subroutine
end program
