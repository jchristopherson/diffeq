# diffeq
A modern Fortran library providing an object-oriented approach to solving and exploring ordinary differential equations.

## Build Status
[![CMake](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/diffeq/workflows/fpm/badge.svg)](https://github.com/jchristopherson/diffeq/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/diffeq/).

## Available Integrators
- Runge-Kutta, 5th Order with 4th-Order Embedded Estimate (Dormand-Prince)
- Tsitouras, 5th Order with 4th-Order Embedded Estimate
- Runge-Kutta, 3rd Order with 2nd-Order Embedded Estimate (Bogacki-Shampine)
- Runge-Kutta, 8th Order with 5th & 3rd-Order Embedded Estimates (Hairer, Nörsett, & Wanner)
- Rosenbrock, 4th Order
- Kennedy-Carpenter ESDIRK, 4th Order with 3rd-Order Embedded Estimate
- Kennedy-Carpenter ESDIRK, 5th Order with 4th-Order Embedded Estimate
- Adams (VODE)
- Backward Differentiation Formula (VODE)

## Getting Started
DIFFEQ solves an initial-value problem of the form

$$
\frac{dy}{dx} = f(x,y), \qquad y(x_0) = y_0.
$$

Define the right-hand side through an `ode_container`, choose an integrator,
and call `solve`. The first column returned by `get_solution` contains the
independent-variable values; the remaining columns contain the solution
components.

```fortran
program first_diffeq
    use iso_fortran_env, only : real64
    use diffeq
    implicit none

    type(ode_container) :: model
    type(runge_kutta_45) :: solver
    real(real64), allocatable :: solution(:,:)

    model%fcn => exponential_rhs
    call solver%solve(model, [0.0d0, 1.0d0], [1.0d0])
    solution = solver%get_solution()

    print *, "y(1) = ", solution(size(solution, 1), 2)

contains
    subroutine exponential_rhs(x, y, dydx, args)
        real(real64), intent(in) :: x
        real(real64), intent(in) :: y(:)
        real(real64), intent(out) :: dydx(:)
        class(*), intent(inout), optional :: args

        dydx(1) = y(1)
    end subroutine exponential_rhs
end program first_diffeq
```

The example solves $y' = y$ with $y(0) = 1$, so the final value is close to
$y(1) = e$. The callback must fill every element of `dydx`; its `args`
argument can be used to pass model parameters to `solve`.

## API Overview
The public API is organized around a model container and interchangeable
integrators:

| API | Purpose |
| --- | --- |
| `ode_container` | Stores the right-hand-side callback and optional Jacobian or mass-matrix callbacks. |
| `model%fcn` | Defines $f(x,y)$ for the problem. This callback is required. |
| `model%jacobian` | Supplies $J = \partial f / \partial y$. If omitted, DIFFEQ estimates it by finite differences where needed. |
| `model%mass_matrix` | Supplies $M(x,y)$ for a system written as $M y' = f(x,y)$. It is supported by the Rosenbrock, Kennedy-Carpenter, and BDF solvers. |
| `integrator%solve(model, x, y0)` | Integrates from `x(1)` to `x(size(x))` using initial state `y0`. |
| `integrator%get_solution()` | Returns an $N \times (n+1)$ array with the independent variable in column 1 and state values in columns 2 through $n+1$. |
| `set_absolute_tolerance` / `set_relative_tolerance` | Control the local error scale, approximately $\mathrm{atol} + \mathrm{rtol}|y|$. |
| `set_maximum_step_size` / `set_minimum_step_size` | Bound adaptive step sizes. |
| `set_allow_overshoot` | Controls whether a final integration step may pass the requested endpoint and be interpolated back. |

The `x` argument may contain only the start and end points, in which case the
solver returns accepted step endpoints. When it contains more than two values,
the solver returns values at those requested points using dense output. All
solver types expose the same `solve` and `get_solution` workflow.

### Choosing an Integrator
- `runge_kutta_23`: inexpensive, lower-order integration for modest accuracy.
- `runge_kutta_45`: general-purpose adaptive integration for non-stiff systems.
- `tsitouras_54`: efficient fifth-order adaptive integration for non-stiff systems; its FSAL tableau uses a fourth-order embedded error estimate.
- `runge_kutta_853`: high-order integration for smooth, non-stiff problems.
- `rosenbrock`: linearly implicit integration for stiff systems and supported
  mass-matrix problems.
- `kennedy_carpenter_4`: fourth-order ESDIRK integration for stiff systems,
    with a third-order embedded error estimate and mass-matrix support.
- `kennedy_carpenter_5`: fifth-order ESDIRK integration for stiff systems,
    with a fourth-order embedded error estimate and mass-matrix support.
- `adams`: variable-order VODE method for smooth, non-stiff systems.
- `bdf`: variable-order VODE method for stiff systems, including supported
    mass-matrix problems.

## Building DIFFEQ
DIFFEQ can be built with either [CMake](https://cmake.org/) or the [Fortran
Package Manager (FPM)](https://github.com/fortran-lang/fpm). Both build systems
require a Fortran 2018 compiler, Git, and BLAS/LAPACK libraries. CMake can
download the LINALG dependency and the test helper automatically when they are
not already installed; FPM resolves the dependencies listed in `fpm.toml`.

### CMake
From the repository root, configure an out-of-source build and compile the
library:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The default CMake build is a static library. To build a shared library instead,
add `-DBUILD_SHARED_LIBS=ON` during configuration. To build the test executable
and register it with CTest, enable `BUILD_TESTING`:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --config Release
ctest --test-dir build --output-on-failure
```

The example programs are enabled separately:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON -DBUILD_DIFFEQ_EXAMPLES=ON
cmake --build build --config Release
```

Install the library, generated Fortran module files, CMake package files, and
pkg-config metadata to a local prefix with:

```sh
cmake --install build --prefix "$PWD/install"
```

On Windows PowerShell, use `cmake --install build --prefix "$PWD/install"`.
For a multi-configuration generator such as Visual Studio, specify the
configuration when installing:

```sh
cmake --install build --config Release --prefix "$PWD/install"
```

To use an installed CMake package from another project, point CMake at the
installation prefix with `-DCMAKE_PREFIX_PATH=/path/to/install` and link the
exported `diffeq::diffeq` target.

### FPM
The repository includes an `fpm.toml` manifest. Build the library with:

```sh
fpm build --profile release
```

Run the test target defined by the manifest with:

```sh
fpm test --profile release
```

FPM places build artifacts under `build/`. The manifest builds DIFFEQ as a
library and disables automatic discovery of executables, examples, and tests;
the declared `diffeq_tests` target is therefore the supported FPM test entry
point.

To use DIFFEQ as a dependency in another FPM project, add it to that project's
`fpm.toml`:

```toml
[dependencies]
diffeq = { git = "https://github.com/jchristopherson/diffeq" }
```

Then import the public module in Fortran with `use diffeq`. The dependency's
LINALG and BLAS/LAPACK requirements are resolved through the consuming
project's FPM and system toolchain configuration.

## Examples
The following example illustrates solving the Van der Pol equation using a 4th-order Rosenbrock solver, but other solvers can be used in an identical manner.  The example also utilizes the [FPLOT](https://github.com/jchristopherson/fplot) library to plot the solution.
```fortran
program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    real(real64), allocatable :: sol(:,:)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd

    ! Define the model
    mdl%fcn => vanderpol

    ! Compute the solution
    call integrator%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
    sol = integrator%get_solution()

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("x")
    call yAxis%set_title("y(x)")
    call y2Axis%set_title("y'(x)")
    call plt%set_use_y2_axis(.true.)
    call lgnd%set_is_visible(.true.)
    
    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_name("y(x)")
    call plt%push(pd1)

    call pd2%define_data(sol(:,1), sol(:,3))
    call pd2%set_draw_against_y2(.true.)
    call pd2%set_name("y'(x)")
    call plt%push(pd2)

    call plt%draw()
end program
```
The routine containing the ODE is located in a different module for this example.  This routine is as follows.
```fortran
pure subroutine vanderpol(x, y, dydx, args)
    ! Arguments
    real(real64), intent(in) :: x, y(:)
    real(real64), intent(out) :: dydx(:)
    class(*), intent(inout), optional :: args

    ! Model Constants
    real(real64), parameter :: mu = 5.0d0

    ! An alternative approach to defining model parameters in the routine is
    ! to use the optional "args" variable.  For example, if mu were to be
    ! passed from the calling routine the args parameter could be used as 
    ! follows:
    !
    ! select type (args)
    ! type is (real(real64))
    !   mu = args
    ! end select
    !
    ! Of course, the call to the solve routine would have to pass this
    ! argument in a manner similar to the following.
    !
    ! mu = 5.0d0
    ! call integrator%solve(...., args = mu)

    ! Equations
    dydx(1) = y(2)
    dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
end subroutine
```
![](images/rosenbrock_example.png?raw=true)



Here's another example comparing the behavior of several integrators for the same Van der Pol problem illustrated in the previous example.  In this example it can be seen that all of the integrators can be utilized in an identical manner.  Additionally, this example illustrates the use of a PI-type controller for step-size control.  Such a controller can be beneficial in the event stability issues are encountered during solution; however, this benefit usually comes with a drawback of decreased efficiency.  For this reason, the default behavior for any of the solvers is to not utilize any PI control; however, it is available if needed.
```fortran
program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    implicit none

    ! Initial Conditions & Time Constraints
    real(real64), parameter :: t(2) = [0.0d0, 5.0d1]
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]

    ! Local Variables
    type(runge_kutta_23) :: integrator_1
    type(runge_kutta_45) :: integrator_2
    type(runge_kutta_853) :: integrator_3
    type(rosenbrock) :: integrator_4
    type(bdf) :: integrator_5
    type(adams) :: integrator_6
    type(kennedy_carpenter_4) :: integrator_7
    type(kennedy_carpenter_5) :: integrator_8
    type(tsitouras_54) :: integrator_9
    type(ode_container) :: mdl
    real(real64), allocatable, dimension(:,:) :: s1, s2, s3, s4, s4a, s5, s6, &
        s7, s8, s9

    ! Define the model
    mdl%fcn => vanderpol

    ! Integrate the model with each integrator
    call integrator_1%solve(mdl, t, ic)
    call integrator_2%solve(mdl, t, ic)
    call integrator_3%solve(mdl, t, ic)
    call integrator_4%solve(mdl, t, ic)
    call integrator_5%solve(mdl, t, ic)
    call integrator_6%solve(mdl, t, ic)
    call integrator_7%solve(mdl, t, ic)
    call integrator_8%solve(mdl, t, ic)
    call integrator_9%solve(mdl, t, ic)

    ! Retrieve the solution from each integrator
    s1 = integrator_1%get_solution()
    s2 = integrator_2%get_solution()
    s3 = integrator_3%get_solution()
    s4 = integrator_4%get_solution()
    s5 = integrator_5%get_solution()
    s6 = integrator_6%get_solution()
    s7 = integrator_7%get_solution()
    s8 = integrator_8%get_solution()
    s9 = integrator_9%get_solution()

    ! Print out the size of each solution
    print "(A, I0, A)", "RUNGE_KUTTA_23: ", size(s1, 1), " Solution Points"
    print "(A, I0, A)", "RUNGE_KUTTA_45: ", size(s2, 1), " Solution Points"
    print "(A, I0, A)", "RUNGE_KUTTA_853: ", size(s3, 1), " Solution Points"
    print "(A, I0, A)", "ROSENBROCK: ", size(s4, 1), " Solution Points"

    ! Now, implement a PI controller and check its effect.  This will likely
    ! increase the number of steps (loss of efficiency), but if there were
    ! any stability issues, stability will likely improve.  Stability is likely
    ! not relevant on this problem, but it's here for illustration purposes.
    call integrator_4%clear_buffer()
    call integrator_4%set_step_size_control_parameter(0.1d0)
    call integrator_4%solve(mdl, t, ic)
    s4a = integrator_4%get_solution()
    print "(A, I0 ,A)", "ROSENBROCK w/ PI Controller: ", size(s4a, 1), " Solution Points"

    ! VODE Integrators
    print "(A, I0, A)", "BDF: ", size(s5, 1), " Solution Points"
    print "(A, I0, A)", "ADAMS: ", size(s6, 1), " Solution Points"

    ! Kennedy-Carpenter Integrators
    print "(A, I0, A)", "KC4: ", size(s7, 1), " Solution Points"
    print "(A, I0, A)", "KC5: ", size(s8, 1), " Solution Points"

    ! Tsitouras Integrators
    print "(A, I0, A)", "TSITOURAS 4/5: ", size(s9, 1), " Solution Points"
end program
```
```txt
RUNGE_KUTTA_23: 2465 Solution Points
RUNGE_KUTTA_45: 583 Solution Points
RUNGE_KUTTA_853: 925 Solution Points
ROSENBROCK: 1187 Solution Points
ROSENBROCK w/ PI Controller: 1187 Solution Points
BDF: 1527 Solution Points
ADAMS: 1865 Solution Points
KC4: 483 Solution Points
KC5: 245 Solution Points
TSITOURAS 4/5: 522 Solution Points
```

## External Libraries
Here is a list of external code libraries utilized by this library.  The CMake build script will include these dependencies automatically; however, it is highly recommended that an optimized BLAS and LAPACK already reside on your system for best performance (used by LINALG for linear algebra calculations).
- [LINALG](https://github.com/jchristopherson/linalg)

## References
1. Butcher, J. C. (2003). Numerical methods for ordinary differential equations. J. Wiley.
2. Shampine, L. F., & Reichelt, M. W., (1997). The MATLAB ODE suite. SIAM Journal on Scientific Computing. 18. 10.1137/S1064827594276424. 
3. J.R. Dormand, P.J. Prince (1980). A family of embedded Runge-Kutta formulae, Journal of Computational and Applied Mathematics, Volume 6, Issue 1, Pages 19-26, ISSN 0377-0427, https://doi.org/10.1016/0771-050X(80)90013-3.
4. P. Bogacki, L.F. Shampine (1989). A 3(2) pair of Runge - Kutta formulas, Applied Mathematics Letters, Volume 2, Issue 4, Pages 321-325, ISSN 0893-9659, https://doi.org/10.1016/0893-9659(89)90079-7.
5. Dormand, J. R. (1996). Numerical methods for differential equations a computational approach. CRC Press. 
6. Kennedy, C. A., Carpenter, M. H. (2016, March). Diagonally Implicit Runge-Kutta Methods for Ordinary Differential Equations. A Review. https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf 
7. Stal, J. (2015). Implementation of Singly Diagonally Implicit Runge-Kutta Methods with Constant Step Sizes. https://core.ac.uk/download/pdf/289938621.pdf 
8. Nayfeh, A. H., & Balachandran, B. (1995). Applied Nonlinear Dynamics: Analytical, Computational, and Experimental Methods. J. Wiley.
