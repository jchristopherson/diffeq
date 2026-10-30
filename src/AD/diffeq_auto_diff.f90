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

module diffeq_auto_diff
    !! First-order forward automatic differentiation for Jacobian expressions.
    !!
    !! Each real_ad stores one value and one directional derivative. Reevaluate
    !! an expression with one independent variable seeded at a time to obtain
    !! Jacobian columns. All overloads are elemental and preserve scalar expansion.
    !! Arithmetic and comparisons accept real_ad, real(real64), and default
    !! integer operands; binary intrinsics accept real_ad and real(real64).
    !!
    !! Comparisons ignore derivatives. ABS uses derivative zero at zero; MIN/MAX
    !! are binary and retain the first argument on ties. SIGN differentiates only
    !! its magnitude argument. These are branch conventions, not derivatives at
    !! discontinuities. Other domain errors and singular derivatives follow the
    !! underlying real intrinsics and the caller's floating-point error settings.
    !!
    !! In the procedures below, a/b (or y/x for ATAN2) are input operands, c is
    !! the result, and dx is propagated by the chain rule. Real and integer
    !! operands are constants with derivative zero.
    use iso_fortran_env
    implicit none
    private
    public :: real_ad
    public :: assignment(=)
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    public :: operator(==), operator(/=), operator(<), operator(<=)
    public :: operator(>), operator(>=)
    public :: abs, sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
    public :: atan2, sinh, cosh, tanh, asinh, acosh, atanh, erf, erfc
    public :: min, max, sign

    type real_ad
        !! A value and its directional derivative. Seed dx with one for an
        !! independent variable and zero for constants.
        real(real64) :: x
            !! Function or independent-variable value.
        real(real64) :: dx
            !! Directional derivative with respect to the current seed.
    end type

    interface assignment(=)
        !! Copy AD data, assign a constant (clearing dx), or extract a real value.
        module procedure :: ad_assign_1
        module procedure :: ad_assign_2
        module procedure :: ad_assign_3
        module procedure :: ad_assign_integer
    end interface

    interface operator(+)
        !! Unary identity and addition, including constant operands in either order.
        module procedure :: ad_add, ad_add_real, real_add_ad
        module procedure :: ad_add_integer, integer_add_ad, ad_positive
    end interface
    interface operator(-)
        !! Unary negation and subtraction, including constant operands.
        module procedure :: ad_subtract, ad_subtract_real, real_subtract_ad
        module procedure :: ad_subtract_integer, integer_subtract_ad, ad_negative
    end interface
    interface operator(*)
        !! Multiplication using the product rule.
        module procedure :: ad_multiply, ad_multiply_real, real_multiply_ad
        module procedure :: ad_multiply_integer, integer_multiply_ad
    end interface
    interface operator(/)
        !! Division using the quotient rule; the denominator must be nonzero.
        module procedure :: ad_divide, ad_divide_real, real_divide_ad
        module procedure :: ad_divide_integer, integer_divide_ad
    end interface
    interface operator(**)
        !! Exponentiation with constant or variable base/exponent.
        !! Variable exponents require a positive base, except that a zero constant
        !! base with a positive exponent has zero derivative. Constant integral
        !! exponents also support negative bases.
        module procedure :: ad_power, ad_power_real, real_power_ad
        module procedure :: ad_power_integer, integer_power_ad
    end interface
    interface operator(==)
        !! Equality of values, independent of derivative seeds.
        module procedure :: ad_equal, ad_equal_real, real_equal_ad
        module procedure :: ad_equal_integer, integer_equal_ad
    end interface
    interface operator(/=)
        !! Inequality of values, independent of derivative seeds.
        module procedure :: ad_not_equal, ad_not_equal_real, real_not_equal_ad
        module procedure :: ad_not_equal_integer, integer_not_equal_ad
    end interface
    interface operator(<)
        !! Strict less-than comparison of values.
        module procedure :: ad_less, ad_less_real, real_less_ad
        module procedure :: ad_less_integer, integer_less_ad
    end interface
    interface operator(<=)
        !! Less-than-or-equal comparison of values.
        module procedure :: ad_less_equal, ad_less_equal_real, real_less_equal_ad
        module procedure :: ad_less_equal_integer, integer_less_equal_ad
    end interface
    interface operator(>)
        !! Strict greater-than comparison of values.
        module procedure :: ad_greater, ad_greater_real, real_greater_ad
        module procedure :: ad_greater_integer, integer_greater_ad
    end interface
    interface operator(>=)
        !! Greater-than-or-equal comparison of values.
        module procedure :: ad_greater_equal, ad_greater_equal_real, real_greater_equal_ad
        module procedure :: ad_greater_equal_integer, integer_greater_equal_ad
    end interface
    interface abs
        !! Absolute value; the derivative at zero is defined to be zero.
        module procedure :: ad_abs
    end interface
    interface sqrt
        !! Square root for nonnegative values; its active derivative is singular at zero.
        module procedure :: ad_sqrt
    end interface
    interface exp
        !! Natural exponential.
        module procedure :: ad_exp
    end interface
    interface log
        !! Natural logarithm for positive values.
        module procedure :: ad_log
    end interface
    interface log10
        !! Base-ten logarithm for positive values.
        module procedure :: ad_log10
    end interface
    interface sin
        !! Sine, with arguments in radians.
        module procedure :: ad_sin
    end interface
    interface cos
        !! Cosine, with arguments in radians.
        module procedure :: ad_cos
    end interface
    interface tan
        !! Tangent, with arguments in radians and away from its poles.
        module procedure :: ad_tan
    end interface
    interface asin
        !! Inverse sine in radians; the active derivative is singular at +/-1.
        module procedure :: ad_asin
    end interface
    interface acos
        !! Inverse cosine in radians; the active derivative is singular at +/-1.
        module procedure :: ad_acos
    end interface
    interface atan
        !! One- or two-argument inverse tangent in radians.
        module procedure :: ad_atan, ad_atan2, ad_atan2_real, real_atan2_ad
    end interface
    interface atan2
        !! Quadrant-aware inverse tangent ATAN2(y,x); not defined at the origin.
        module procedure :: ad_atan2, ad_atan2_real, real_atan2_ad
    end interface
    interface sinh
        !! Hyperbolic sine.
        module procedure :: ad_sinh
    end interface
    interface cosh
        !! Hyperbolic cosine.
        module procedure :: ad_cosh
    end interface
    interface tanh
        !! Hyperbolic tangent.
        module procedure :: ad_tanh
    end interface
    interface asinh
        !! Inverse hyperbolic sine.
        module procedure :: ad_asinh
    end interface
    interface acosh
        !! Inverse hyperbolic cosine for values >= 1; active derivative singular at 1.
        module procedure :: ad_acosh
    end interface
    interface atanh
        !! Inverse hyperbolic tangent for values strictly between -1 and 1.
        module procedure :: ad_atanh
    end interface
    interface erf
        !! Error function.
        module procedure :: ad_erf
    end interface
    interface erfc
        !! Complementary error function.
        module procedure :: ad_erfc
    end interface
    interface min
        !! Binary minimum, carrying the selected derivative; first argument wins ties.
        module procedure :: ad_min, ad_min_real, real_min_ad
    end interface
    interface max
        !! Binary maximum, carrying the selected derivative; first argument wins ties.
        module procedure :: ad_max, ad_max_real, real_max_ad
    end interface
    interface sign
        !! Magnitude of a with the sign of b; ignores the derivative of b.
        module procedure :: ad_sign, ad_sign_real, real_sign_ad
    end interface

contains
! ******************************************************************************
! ASSIGNMENT OPERATORS
! ------------------------------------------------------------------------------
elemental subroutine ad_assign_1(lhs, rhs)
    !! Copy the value and derivative of rhs to lhs.
    type(real_ad), intent(out) :: lhs
    type(real_ad), intent(in) :: rhs

    lhs%x = rhs%x
    lhs%dx = rhs%dx
end subroutine

elemental subroutine ad_assign_2(lhs, rhs)
    !! Assign a real64 constant rhs to lhs and clear its derivative.
    type(real_ad), intent(out) :: lhs
    real(real64), intent(in) :: rhs
    lhs%x = rhs
    lhs%dx = 0.0d0
end subroutine

elemental subroutine ad_assign_3(lhs, rhs)
    !! Extract the value of rhs into real64 lhs, discarding the derivative.
    real(real64), intent(out) :: lhs
    type(real_ad), intent(in) :: rhs
    lhs = rhs%x
end subroutine

elemental subroutine ad_assign_integer(lhs, rhs)
    !! Assign a default-integer constant rhs to lhs and clear its derivative.
    type(real_ad), intent(out) :: lhs
    integer, intent(in) :: rhs
    lhs%x = real(rhs, real64)
    lhs%dx = 0.0d0
end subroutine

! ******************************************************************************
! MATHEMATICAL OPERATORS
! ------------------------------------------------------------------------------
elemental function ad_positive(a) result(c)
    !! Unary plus preserves both components.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = a
end function

elemental function ad_negative(a) result(c)
    !! Unary minus negates both components.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(-a%x, -a%dx)
end function

elemental function ad_add(a, b) result(c)
    !! Add values and directional derivatives.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c = real_ad(a%x + b%x, a%dx + b%dx)
end function

elemental function ad_add_real(a, b) result(c)
    !! Add a real64 constant without changing the derivative.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(a%x + b, a%dx)
end function

elemental function real_add_ad(a, b) result(c)
    !! Add an AD operand to a real64 constant.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = b + a
end function

elemental function ad_add_integer(a, b) result(c)
    !! Add a default-integer constant using real64 arithmetic.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    type(real_ad) :: c
    c = a + real(b, real64)
end function

elemental function integer_add_ad(a, b) result(c)
    !! Add an AD operand to a default-integer constant.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real(a, real64) + b
end function

elemental function ad_subtract(a, b) result(c)
    !! Subtract values and directional derivatives.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c = real_ad(a%x - b%x, a%dx - b%dx)
end function

elemental function ad_subtract_real(a, b) result(c)
    !! Subtract a real64 constant without changing the derivative.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(a%x - b, a%dx)
end function

elemental function real_subtract_ad(a, b) result(c)
    !! Subtract an AD operand from a real64 constant, negating its derivative.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(a - b%x, -b%dx)
end function

elemental function ad_subtract_integer(a, b) result(c)
    !! Subtract a default-integer constant using real64 arithmetic.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    type(real_ad) :: c
    c = a - real(b, real64)
end function

elemental function integer_subtract_ad(a, b) result(c)
    !! Subtract an AD operand from a default-integer constant.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real(a, real64) - b
end function

elemental function ad_multiply(a, b) result(c)
    !! Apply the product rule: d(a*b) = da*b + a*db.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c = real_ad(a%x * b%x, a%dx * b%x + a%x * b%dx)
end function

elemental function ad_multiply_real(a, b) result(c)
    !! Scale the value and derivative by a real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(a%x * b, a%dx * b)
end function

elemental function real_multiply_ad(a, b) result(c)
    !! Multiply a real64 constant by an AD operand.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = b * a
end function

elemental function ad_multiply_integer(a, b) result(c)
    !! Scale an AD operand by a default-integer constant.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    type(real_ad) :: c
    c = a * real(b, real64)
end function

elemental function integer_multiply_ad(a, b) result(c)
    !! Multiply a default-integer constant by an AD operand.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real(a, real64) * b
end function

elemental function ad_divide(a, b) result(c)
    !! Apply the quotient rule as d(a/b) = (da - (a/b)*db)/b.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c%x = a%x / b%x
    c%dx = (a%dx - c%x * b%dx) / b%x
end function

elemental function ad_divide_real(a, b) result(c)
    !! Divide the value and derivative by a nonzero real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(a%x / b, a%dx / b)
end function

elemental function real_divide_ad(a, b) result(c)
    !! Divide a real64 constant by AD b, with derivative -a*db/b**2.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c%x = a / b%x
    c%dx = -c%x * (b%dx / b%x)
end function

elemental function ad_divide_integer(a, b) result(c)
    !! Divide an AD operand by a nonzero default-integer constant.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    type(real_ad) :: c
    c = a / real(b, real64)
end function

elemental function integer_divide_ad(a, b) result(c)
    !! Divide a default-integer constant by an AD operand.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real(a, real64) / b
end function

elemental function ad_power(a, b) result(c)
    !! Apply d(a**b) = b*a**(b-1)*da + a**b*log(a)*db.
    !! Omit inactive terms so constant powers work at zero and negative bases.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c = a ** b%x
    ! A constant exponent must not introduce log(a), which would exclude
    ! negative bases even for integral powers.
    if (b%dx /= 0.0d0) then
        if (a%x /= 0.0d0 .or. b%x <= 0.0d0) then
            c%dx = c%dx + c%x * log(a%x) * b%dx
        end if
    end if
end function

elemental function ad_power_real(a, b) result(c)
    !! Raise AD a to constant b with derivative b*a**(b-1)*da.
    !! Exponents zero and one are handled without singular intermediate powers.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c%x = a%x ** b
    if (b == 0.0d0 .or. a%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (b == 1.0d0) then
        c%dx = a%dx
    else
        c%dx = b * (a%x ** (b - 1.0d0)) * a%dx
    end if
end function

elemental function real_power_ad(a, b) result(c)
    !! Raise constant a to AD b with derivative a**b*log(a)*db.
    !! A zero base and positive exponent give zero for both components.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c%x = a ** b%x
    if (b%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (a == 0.0d0 .and. b%x > 0.0d0) then
        c%dx = 0.0d0
    else
        c%dx = c%x * log(a) * b%dx
    end if
end function

elemental function ad_power_integer(a, b) result(c)
    !! Evaluate an integer power, preserving negative-base and zero-base behavior.
    !! Negative exponents avoid b-1, which could overflow the integer range.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    type(real_ad) :: c
    c%x = a%x ** b
    if (b == 0 .or. a%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (b > 0) then
        c%dx = real(b, real64) * (a%x ** (b - 1)) * a%dx
    else
        c%dx = real(b, real64) * (c%x / a%x) * a%dx
    end if
end function

elemental function integer_power_ad(a, b) result(c)
    !! Raise a default-integer constant to an AD exponent using real64 arithmetic.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real(a, real64) ** b
end function

! ******************************************************************************
! COMPARISON OPERATORS
! ------------------------------------------------------------------------------
elemental logical function ad_equal(a, b) result(c)
    !! Compare only the values of two AD operands for equality.
    type(real_ad), intent(in) :: a, b
    c = a%x == b%x
end function

elemental logical function ad_equal_real(a, b) result(c)
    !! Compare an AD value with a real64 constant for equality.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x == b
end function

elemental logical function real_equal_ad(a, b) result(c)
    !! Compare a real64 constant with an AD value for equality.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a == b%x
end function

elemental logical function ad_equal_integer(a, b) result(c)
    !! Compare an AD value with a default integer for equality.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x == real(b, real64)
end function

elemental logical function integer_equal_ad(a, b) result(c)
    !! Compare a default integer with an AD value for equality.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) == b%x
end function

elemental logical function ad_not_equal(a, b) result(c)
    !! Compare only the values of two AD operands for inequality.
    type(real_ad), intent(in) :: a, b
    c = a%x /= b%x
end function

elemental logical function ad_not_equal_real(a, b) result(c)
    !! Compare an AD value with a real64 constant for inequality.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x /= b
end function

elemental logical function real_not_equal_ad(a, b) result(c)
    !! Compare a real64 constant with an AD value for inequality.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a /= b%x
end function

elemental logical function ad_not_equal_integer(a, b) result(c)
    !! Compare an AD value with a default integer for inequality.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x /= real(b, real64)
end function

elemental logical function integer_not_equal_ad(a, b) result(c)
    !! Compare a default integer with an AD value for inequality.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) /= b%x
end function

elemental logical function ad_less(a, b) result(c)
    !! Test whether AD a's value is strictly less than AD b's.
    type(real_ad), intent(in) :: a, b
    c = a%x < b%x
end function

elemental logical function ad_less_real(a, b) result(c)
    !! Test whether an AD value is less than a real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x < b
end function

elemental logical function real_less_ad(a, b) result(c)
    !! Test whether a real64 constant is less than an AD value.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a < b%x
end function

elemental logical function ad_less_integer(a, b) result(c)
    !! Test whether an AD value is less than a default integer.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x < real(b, real64)
end function

elemental logical function integer_less_ad(a, b) result(c)
    !! Test whether a default integer is less than an AD value.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) < b%x
end function

elemental logical function ad_less_equal(a, b) result(c)
    !! Test whether AD a's value is less than or equal to AD b's.
    type(real_ad), intent(in) :: a, b
    c = a%x <= b%x
end function

elemental logical function ad_less_equal_real(a, b) result(c)
    !! Test whether an AD value is at most a real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x <= b
end function

elemental logical function real_less_equal_ad(a, b) result(c)
    !! Test whether a real64 constant is at most an AD value.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a <= b%x
end function

elemental logical function ad_less_equal_integer(a, b) result(c)
    !! Test whether an AD value is at most a default integer.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x <= real(b, real64)
end function

elemental logical function integer_less_equal_ad(a, b) result(c)
    !! Test whether a default integer is at most an AD value.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) <= b%x
end function

elemental logical function ad_greater(a, b) result(c)
    !! Test whether AD a's value is strictly greater than AD b's.
    type(real_ad), intent(in) :: a, b
    c = a%x > b%x
end function

elemental logical function ad_greater_real(a, b) result(c)
    !! Test whether an AD value exceeds a real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x > b
end function

elemental logical function real_greater_ad(a, b) result(c)
    !! Test whether a real64 constant exceeds an AD value.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a > b%x
end function

elemental logical function ad_greater_integer(a, b) result(c)
    !! Test whether an AD value exceeds a default integer.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x > real(b, real64)
end function

elemental logical function integer_greater_ad(a, b) result(c)
    !! Test whether a default integer exceeds an AD value.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) > b%x
end function

elemental logical function ad_greater_equal(a, b) result(c)
    !! Test whether AD a's value is greater than or equal to AD b's.
    type(real_ad), intent(in) :: a, b
    c = a%x >= b%x
end function

elemental logical function ad_greater_equal_real(a, b) result(c)
    !! Test whether an AD value is at least a real64 constant.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    c = a%x >= b
end function

elemental logical function real_greater_equal_ad(a, b) result(c)
    !! Test whether a real64 constant is at least an AD value.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    c = a >= b%x
end function

elemental logical function ad_greater_equal_integer(a, b) result(c)
    !! Test whether an AD value is at least a default integer.
    type(real_ad), intent(in) :: a
    integer, intent(in) :: b
    c = a%x >= real(b, real64)
end function

elemental logical function integer_greater_equal_ad(a, b) result(c)
    !! Test whether a default integer is at least an AD value.
    integer, intent(in) :: a
    type(real_ad), intent(in) :: b
    c = real(a, real64) >= b%x
end function

! ******************************************************************************
! INTRINSIC FUNCTIONS
! ------------------------------------------------------------------------------
elemental function ad_abs(a) result(c)
    !! Propagate sign(a)*da, choosing derivative zero at a = 0.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c%x = abs(a%x)
    if (a%x == 0.0d0) then
        c%dx = 0.0d0
    else
        c%dx = sign(1.0d0, a%x) * a%dx
    end if
end function

elemental function ad_sqrt(a) result(c)
    !! Propagate da/(2*sqrt(a)); inactive seeds remain zero at a = 0.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c%x = sqrt(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = 0.5d0 * a%dx / c%x
end function

elemental function ad_exp(a) result(c)
    !! Propagate exp(a)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c%x = exp(a%x)
    c%dx = c%x * a%dx
end function

elemental function ad_log(a) result(c)
    !! Propagate da/a for the natural logarithm.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(log(a%x), a%dx / a%x)
end function

elemental function ad_log10(a) result(c)
    !! Propagate da/(a*log(10)) for the base-ten logarithm.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(log10(a%x), (a%dx / a%x) / log(10.0d0))
end function

elemental function ad_sin(a) result(c)
    !! Propagate cos(a)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(sin(a%x), cos(a%x) * a%dx)
end function

elemental function ad_cos(a) result(c)
    !! Propagate -sin(a)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(cos(a%x), -sin(a%x) * a%dx)
end function

elemental function ad_tan(a) result(c)
    !! Propagate da/cos(a)**2.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(tan(a%x), a%dx / cos(a%x)**2)
end function

elemental function ad_asin(a) result(c)
    !! Propagate da/sqrt(1-a**2), omitting the singular term for inactive seeds.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c%x = asin(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = a%dx / sqrt((1.0d0 - a%x) * (1.0d0 + a%x))
end function

elemental function ad_acos(a) result(c)
    !! Propagate -da/sqrt(1-a**2) using the inverse-sine derivative.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = asin(a)
    c%x = acos(a%x)
    c%dx = -c%dx
end function

elemental function ad_atan(a) result(c)
    !! Propagate da/(1+a**2), avoiding overflow for large a.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    real(real64) :: inverse_x
    c%x = atan(a%x)
    if (abs(a%x) <= 1.0d0) then
        c%dx = a%dx / (1.0d0 + a%x * a%x)
    else
        inverse_x = 1.0d0 / a%x
        c%dx = ((a%dx / a%x) / a%x) / (1.0d0 + inverse_x * inverse_x)
    end if
end function

elemental function ad_atan2(y, x) result(c)
    !! Propagate (x*dy-y*dx)/(x**2+y**2) using normalized coordinates.
    type(real_ad), intent(in) :: y, x
    type(real_ad) :: c
    real(real64) :: scale, xscaled, yscaled
    scale = max(abs(x%x), abs(y%x))
    xscaled = x%x / scale
    yscaled = y%x / scale
    c%x = atan2(y%x, x%x)
    c%dx = (xscaled * (y%dx / scale) - yscaled * (x%dx / scale)) / &
        (xscaled * xscaled + yscaled * yscaled)
end function

elemental function ad_atan2_real(y, x) result(c)
    !! Evaluate ATAN2 with AD ordinate y and constant abscissa x.
    type(real_ad), intent(in) :: y
    real(real64), intent(in) :: x
    type(real_ad) :: c
    c = atan2(y, real_ad(x, 0.0d0))
end function

elemental function real_atan2_ad(y, x) result(c)
    !! Evaluate ATAN2 with constant ordinate y and AD abscissa x.
    real(real64), intent(in) :: y
    type(real_ad), intent(in) :: x
    type(real_ad) :: c
    c = atan2(real_ad(y, 0.0d0), x)
end function

elemental function ad_sinh(a) result(c)
    !! Propagate cosh(a)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(sinh(a%x), cosh(a%x) * a%dx)
end function

elemental function ad_cosh(a) result(c)
    !! Propagate sinh(a)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(cosh(a%x), sinh(a%x) * a%dx)
end function

elemental function ad_tanh(a) result(c)
    !! Propagate da/cosh(a)**2 using exp(-abs(a)) to avoid hyperbolic overflow.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    real(real64) :: t, sech
    t = exp(-abs(a%x))
    sech = 2.0d0 * t / (1.0d0 + t*t)
    c%x = tanh(a%x)
    c%dx = (a%dx * sech) * sech
end function

elemental function ad_asinh(a) result(c)
    !! Propagate da/sqrt(a**2+1), using HYPOT to avoid overflow.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(asinh(a%x), a%dx / hypot(a%x, 1.0d0))
end function

elemental function ad_acosh(a) result(c)
    !! Propagate da/(sqrt(a-1)*sqrt(a+1)); inactive seeds stay zero at a = 1.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c%x = acosh(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = (a%dx / sqrt(a%x - 1.0d0)) / sqrt(a%x + 1.0d0)
end function

elemental function ad_atanh(a) result(c)
    !! Propagate da/(1-a**2) using factored divisions.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(atanh(a%x), (a%dx / (1.0d0 - a%x)) / (1.0d0 + a%x))
end function

elemental function ad_erf(a) result(c)
    !! Propagate 2/sqrt(pi)*exp(-a**2)*da.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = real_ad(erf(a%x), 2.0d0 / sqrt(acos(-1.0d0)) * exp(-a%x**2) * a%dx)
end function

elemental function ad_erfc(a) result(c)
    !! Negate the error-function derivative for the complementary error function.
    type(real_ad), intent(in) :: a
    type(real_ad) :: c
    c = erf(a)
    c%x = erfc(a%x)
    c%dx = -c%dx
end function

elemental function ad_min(a, b) result(c)
    !! Return the lower-valued operand with its derivative, choosing a on a tie.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    if (a%x <= b%x) then
        c = a
    else
        c = b
    end if
end function

elemental function ad_min_real(a, b) result(c)
    !! Select the minimum of AD a and constant b; selecting b clears the derivative.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = min(a, real_ad(b, 0.0d0))
end function

elemental function real_min_ad(a, b) result(c)
    !! Select the minimum of constant a and AD b; selecting a clears the derivative.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = min(real_ad(a, 0.0d0), b)
end function

elemental function ad_max(a, b) result(c)
    !! Return the higher-valued operand with its derivative, choosing a on a tie.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    if (a%x >= b%x) then
        c = a
    else
        c = b
    end if
end function

elemental function ad_max_real(a, b) result(c)
    !! Select the maximum of AD a and constant b; selecting b clears the derivative.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = max(a, real_ad(b, 0.0d0))
end function

elemental function real_max_ad(a, b) result(c)
    !! Select the maximum of constant a and AD b; selecting a clears the derivative.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = max(real_ad(a, 0.0d0), b)
end function

elemental function ad_sign(a, b) result(c)
    !! Propagate sign(1,b)*d(abs(a)), treating the sign selection as constant.
    type(real_ad), intent(in) :: a, b
    type(real_ad) :: c
    c = abs(a)
    c%x = sign(a%x, b%x)
    c%dx = sign(1.0d0, b%x) * c%dx
end function

elemental function ad_sign_real(a, b) result(c)
    !! Apply a constant sign to AD a and its absolute-value derivative.
    type(real_ad), intent(in) :: a
    real(real64), intent(in) :: b
    type(real_ad) :: c
    c = sign(a, real_ad(b, 0.0d0))
end function

elemental function real_sign_ad(a, b) result(c)
    !! Apply AD b's sign to a constant magnitude, returning derivative zero.
    real(real64), intent(in) :: a
    type(real_ad), intent(in) :: b
    type(real_ad) :: c
    c = real_ad(sign(a, b%x), 0.0d0)
end function

end module