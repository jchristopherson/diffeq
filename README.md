# diffeq
A modern Fortran library providing an object-oriented approach to solving ordinary differential equations and related analyses.

## Work-In-Progress
This library is currently a work-in-progress.  The API is definitely not in a stable condition.

## Build Status
[![CMake](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/diffeq/).

## Available ODE Solvers
### Fixed Step
- 4th Order Runge-Kutta
- Adams-Bashforth-Moulton

### Variable Step
- Dormand-Prince Runge-Kutta 5th/4th Order
- Bogacki-Shampine Runge-Kutta 3rd/2nd Order

## Additional Operations
- FFT-Based Frequency Response Function Estimation
- Frequency-Sweep-Based Frequency Response Function Estimation

## References
1. Butcher, J. C. (2003). Numerical methods for ordinary differential equations. J. Wiley.
2. Shampine, L. F., & Reichelt, M. W., (1997). The MATLAB ODE suite. SIAM Journal on Scientific Computing. 18. 10.1137/S1064827594276424. 
3. J.R. Dormand, P.J. Prince (1980). A family of embedded Runge-Kutta formulae, Journal of Computational and Applied Mathematics, Volume 6, Issue 1, Pages 19-26, ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90013-3.
4. P. Bogacki, L.F. Shampine (1989). A 3(2) pair of Runge - Kutta formulas, Applied Mathematics Letters, Volume 2, Issue 4, Pages 321-325, ISSN 0893-9659, https://doi.org/10.1016/0893-9659(89)90079-7.