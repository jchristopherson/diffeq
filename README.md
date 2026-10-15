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
- Backward Differentiation Formula (BDF), Variable Order (1-5)
- Adams-Bashforth-Moulton Predictor-Corrector (PECE), Variable Order (1-12)

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

Not every problem arrives in the explicit form shown above. Many mechanical
and electrical systems are more naturally written with a mass matrix,

$$
M(x,y) \frac{dy}{dx} = f(x,y),
$$

where $M$ couples the derivatives to one another. Supplying a routine for $M$
through `model%mass_matrix` is supported by the `rosenbrock`,
Kennedy-Carpenter, `bdf`, and `adams` integrators. If the matrix is constant,
calling `model%set_is_mass_matrix_dependent(.false.)` avoids recomputing it at
every step.

So long as $M$ is nonsingular the problem remains an ordinary differential
equation. If $M$ is singular, however, one or more of its rows carry no
derivative at all and instead impose an algebraic relationship among the
states. The system is then a differential-algebraic equation (DAE) rather
than an ODE. The `rosenbrock` integrator forms and factors
$\frac{1}{\gamma h} M - J$ by means of a QR factorization with column
pivoting, so it is able to accommodate a singular $M$ and solve index-1
DAEs. The Kennedy-Carpenter ESDIRK integrators use implicit stage solves
and a consistent-derivative calculation, and `bdf` factors
$\frac{\alpha_0}{h} M - J$ in the same spirit. Consequently, all of these
integrators can accommodate a singular $M$ and solve index-1 DAEs. The
[Differential-Algebraic Equations](#differential-algebraic-equations) example
works such a problem. Two points are worth noting when doing so. The initial
conditions handed to `solve` must satisfy the algebraic constraints, and a
problem of index greater than one must first be reduced to index 1, typically
by differentiating its constraint equations with respect to the independent
variable.

For autonomous holonomic mechanical constraints, `ggl_container` provides
a [stabilized index-1 GGL variant](#stabilized-index-1-ggl) that assembles
the reduced equations and damps position and velocity constraint violations.

The `adams` integrator is the exception. It integrates $y'$ directly, so it
must resolve the mass matrix rather than merely carry it along as a factor.
It therefore requires $M$ to be nonsingular and reports
`DIFFEQ_SINGULAR_MATRIX_ERROR` if it is not. A DAE must be given to one of
the other integrators.

## API Overview
The public API is organized around a model container and interchangeable
integrators:

| API | Purpose |
| --- | --- |
| `ode_container` | Stores the right-hand-side callback, optional Jacobian, and optional mass-matrix callbacks. |
| `ggl_container` | Adapts holonomically constrained mechanics to a stabilized index-1 DAE with position and velocity residual damping. Its force and physical-mass callbacks have specialized shapes documented below. |
| `model%fcn` | Defines $f(x,y)$ for the problem. This callback is required. |
| `model%jacobian` | Supplies $J = \partial f / \partial y$. If omitted, DIFFEQ estimates it by finite differences where needed. |
| `model%mass_matrix` | Supplies $M(x,y)$ for a system written as $M y' = f(x,y)$. It is supported by the Rosenbrock, Kennedy-Carpenter, BDF, and Adams solvers. Only Adams requires $M$ to be nonsingular. |
| `integrator%solve(model, x, y0)` | Integrates from `x(1)` to `x(size(x))` using initial state `y0`. |
| `integrator%get_solution()` | Returns an $N \times (n+1)$ array with the independent variable in column 1 and state values in columns 2 through $n+1$. |
| `integrator%get_ode_results()` | Returns an `ode_results` container that provides the solution and statistics on integrator performance, including step counts and function, Jacobian, and matrix-factorization evaluations. |
| `set_absolute_tolerance` / `set_relative_tolerance` | Control the local error scale, approximately $\mathrm{atol} + \mathrm{rtol}|y|$. |
| `set_step_size_control_parameter` | Sets the PI controller parameter used by the Kennedy-Carpenter solvers. Rosenbrock and the multi-step solvers use their own step-size estimators. |
| `set_maximum_order` / `get_maximum_order` | Bound the order of a multi-step solver. A request above what the method supports is clamped rather than rejected. |
| `set_newton_step_limit` / `set_newton_tolerance` | Govern the Newton iteration used by the implicit multi-step solvers, such as `bdf`. |
| `set_maximum_step_size` / `set_minimum_step_size` | Bound adaptive step sizes. |
| `set_allow_overshoot` | Controls whether a final integration step may pass the requested endpoint and be interpolated back. |

Tolerance setters accept either a scalar applied to every state variable or
an array containing one value per state variable:

```fortran
call integrator%set_absolute_tolerance(1.0d-8)
call integrator%set_relative_tolerance([1.0d-6, 1.0d-9])
```

An array must have the same size as the initial state passed to `solve`.
Calling the same setter later with a scalar restores the blanket tolerance.
Per-state tolerances apply to every adaptive solver, including Rosenbrock,
Kennedy-Carpenter, BDF, and Adams when they select their next step or order.

The `x` argument may contain only the start and end points, in which case the
solver returns accepted step endpoints. When it contains more than two values,
the solver returns values at those requested points using dense output. All
solver types expose the same `solve` and `get_solution` workflow.

### Choosing an Integrator
- `runge_kutta_23`
    - Capabilities: inexpensive, lower-order adaptive integration for modest
        accuracy and non-stiff systems.
    - Limitations: explicit method; not intended for stiff systems or
        mass-matrix DAEs.
- `runge_kutta_45`
    - Capabilities: general-purpose adaptive integration for non-stiff systems;
        uses an embedded error estimate and FSAL stage reuse.
    - Limitations: explicit method; not intended for stiff systems or
        mass-matrix DAEs.
- `tsitouras_54`
    - Capabilities: efficient fifth-order adaptive integration for non-stiff
        systems; its FSAL tableau uses a fourth-order embedded error estimate.
    - Limitations: explicit method; not intended for stiff systems or
        mass-matrix DAEs.
- `runge_kutta_853`
    - Capabilities: high-order adaptive integration for smooth, non-stiff
        problems, with fifth- and third-order embedded estimates.
    - Limitations: explicit method; not intended for stiff systems or
        mass-matrix DAEs.
- `rosenbrock`
    - Capabilities: linearly implicit integration for stiff systems and supported
        mass-matrix problems. It accommodates singular mass matrices by factoring
        $\frac{1}{\gamma h} M - J$ rather than inverting $M$ on its own, and can
        solve index-1 DAEs when the initial conditions satisfy the algebraic
        constraints.
    - Limitations: uses a Rosenbrock-specific step-size estimator rather than
        the inherited PI controller; higher-index DAEs must first be reduced to
        index 1.
- `kennedy_carpenter_4`
    - Capabilities: fourth-order ESDIRK integration for stiff systems, with a
        third-order embedded error estimate. Supports nonsingular and singular
        mass matrices through diagonally implicit stage solves and the configurable
        inherited PI step-size controller.
    - Limitations: singular-mass problems require consistent initial conditions
        and must be index 1.
- `kennedy_carpenter_5`
    - Capabilities: fifth-order ESDIRK integration for stiff systems, with a
        fourth-order embedded error estimate. Supports nonsingular and singular
        mass matrices through diagonally implicit stage solves and the configurable
        inherited PI step-size controller.
    - Limitations: singular-mass problems require consistent initial conditions
        and must be index 1.
- `bdf`
    - Capabilities: variable step-size, variable order (1-5) backward
        differentiation formulae for stiff systems and index-1 DAEs. The
        coefficients are formed from the actual spacing of the stored solution
        points, so a change of step size is handled exactly. Accommodates
        singular mass matrices by factoring $\frac{\alpha_0}{h} M - J$ rather
        than inverting $M$. Being a multi-step method, it reaches a given order
        on far fewer derivative evaluations than a Runge-Kutta method of the
        same order, and its dense output costs no extra evaluations at all.
    - Limitations: the order is capped at five, as BDF formulae lose zero
        stability beyond six. Each step requires a Jacobian and the solution of
        a linear system. The method restarts at first order, so problems with
        frequent discontinuities forfeit much of its advantage. Singular-mass
        problems require consistent initial conditions and must be index 1.
- `adams`
    - Capabilities: variable step-size, variable order (1-12)
        Adams-Bashforth-Moulton predictor-corrector in PECE mode, for non-stiff
        systems. Each step costs exactly two derivative evaluations and needs
        neither a Jacobian nor a linear solve, which makes it the most
        economical choice when the problem is not stiff and the tolerances are
        tight. Its high order ceiling pays off there: on a smooth problem at a
        tolerance of $10^{-12}$ it reached order nine and used less than half
        the steps the same method restricted to order five required.
    - Limitations: the corrector is applied once rather than iterated, so it is
        contractive only while $hL < 1$ and is unsuited to stiff problems; it
        will still hold its accuracy on one, but only by taking orders of
        magnitude more steps than an implicit method needs. A mass matrix must
        be nonsingular, so DAEs are out of reach. The method restarts at first
        order.

### Step Size and Order Control in the Multi-Step Solvers
`bdf` and `adams` share a single controller, so what follows applies equally to both; only the error constants belonging to each family of formulae differ.

The predictor and the corrector are of the same order, so their difference is proportional to the local truncation error of the step.  Scaling that difference by the tolerances and taking a root-mean-square norm gives a normalized error $e$, and the step is accepted when $e \le 1$.

The order is reconsidered after every accepted step.  Divided differences of the stored solution supply the error the method would have committed at orders $k-1$, $k$, and $k+1$, and each is converted into the step-size ratio it would permit,

$$
\rho_j = \left( \frac{1}{e_j} \right)^{1/(j+1)}.
$$

Whichever order offers the largest $\rho$ is adopted.  Measuring all three candidates the same way matters here: a comparison drawn between differently scaled estimates would favor one order for no better reason than its constant.  An increase in order is withheld until the current order has survived $k+1$ consecutive steps, which keeps the method from chasing transients, and the order can never exceed the history on hand.  Both solvers therefore begin at first order, which is also why neither requires a starting procedure or a consistent initial derivative.

The step size for the next step is the $\rho$ belonging to the chosen order, multiplied by the safety factor.  Growth is capped at tenfold and reduction at one tenth, and a ratio falling between one and two is treated as no change at all: a multi-step formula interpolates its own history, so its stability suffers when the spacing of that history is disturbed without something to gain.  A step failing the error test is retried at between one half and one tenth of its size, a corrector that fails outright quarters it, and three consecutive failures of either kind return the order to one.

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
### The Van der Pol Equation
The Van der Pol equation describes a self-sustaining oscillator, and is written as the second-order equation

$$
\frac{d^2 y}{dx^2} - \mu \left( 1 - y^2 \right) \frac{dy}{dx} + y = 0.
$$

The quantity $\mu \left( 1 - y^2 \right)$ acts as a damping coefficient whose sign depends upon the amplitude of the solution.  When $\left| y \right| > 1$ the quantity is negative and energy is removed from the system, and when $\left| y \right| < 1$ it is positive and energy is returned.  The consequence is that every non-trivial solution is drawn onto the same closed orbit, or limit cycle, regardless of the conditions from which it started.  The parameter $\mu$ governs how abruptly that exchange of energy takes place.  As $\mu$ increases the solution lingers longer on the slow portions of the cycle and transitions between them ever more rapidly, and the problem becomes increasingly stiff.

The integrators in this library operate on systems of first-order equations, so the equation is recast as a pair of first-order equations by letting $y_1 = y$ and $y_2 = dy/dx$.

$$
\dot{y}_1 = y_2, \qquad \dot{y}_2 = \mu \left( 1 - y_1^2 \right) y_2 - y_1.
$$

The following example solves this system for $\mu = 5$ over $0 \le x \le 50$ with $y_1(0) = 2$ and $y_2(0) = 0$ using a 4th-order Rosenbrock solver, but other solvers can be used in an identical manner.  The example also utilizes the [FPLOT](https://github.com/jchristopherson/fplot) library to plot the solution.
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

The limit cycle is readily apparent in the results.  The solution $y(x)$ remains near $\pm 2$ while drifting slowly, and then reverses sign over a comparatively short span of the independent variable.  Each reversal appears in the derivative $y'(x)$ as a sharp spike whose magnitude is several times the amplitude of $y(x)$ itself.  This disparity between the slow and fast portions of the solution is precisely the behavior that an adaptive step-size algorithm is meant to accommodate, as the integrator must take small steps through each transition while being free to take much larger steps elsewhere.  It is also what makes the equation a useful problem against which to compare integrators, which is the subject of the next example.

### Comparing Integrators
Here's another example comparing the behavior of several integrators on the same Van der Pol problem illustrated in the previous example.  Each integrator is used in the same manner; only its declared type changes.  Calling `get_ode_results()` retrieves an `ode_results` container holding both the solution and an `ode_results_statistics` member named `statistics`.  This provides a more informative comparison than the number of stored solution points alone:

- `statistics%successful_steps`: the number of successful integration steps.
- `statistics%failed_steps`: the number of unsuccessful step attempts.
- `statistics%jacobian_evaluations`: the number of Jacobian matrix evaluations, relevant to methods that use a Jacobian.
- `statistics%function_evaluations`: the number of evaluations of the ODE model function.
- `statistics%matrix_factorizations`: the number of matrix factorizations, relevant to methods that solve linear systems.

The example reports solution points, step statistics, and function evaluations for every method, along with Jacobian evaluations and matrix factorizations for the implicit methods.  It also compares `runge_kutta_45` with and without a PI-type step-size controller.  The complete source is available in [comparison_example.f90](examples/comparison_example.f90).
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
    type(kennedy_carpenter_4) :: integrator_5
    type(kennedy_carpenter_5) :: integrator_6
    type(tsitouras_54) :: integrator_7
    type(bdf) :: integrator_8
    type(adams) :: integrator_9
    type(ode_container) :: mdl
    type(ode_results) :: s1, s2, s2a, s3, s4, s5, s6, s7, s8, s9

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
    s1 = integrator_1%get_ode_results()
    s2 = integrator_2%get_ode_results()
    s3 = integrator_3%get_ode_results()
    s4 = integrator_4%get_ode_results()
    s5 = integrator_5%get_ode_results()
    s6 = integrator_6%get_ode_results()
    s7 = integrator_7%get_ode_results()
    s8 = integrator_8%get_ode_results()
    s9 = integrator_9%get_ode_results()

    ! Print out the size of each solution
    print "(A)", "RUNGE_KUTTA_23: "
    print "(A, I0, A)", achar(9), size(s1, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s1%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s1%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s1%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "RUNGE_KUTTA_45: "
    print "(A, I0, A)", achar(9), size(s2, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s2%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s2%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s2%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "RUNGE_KUTTA_853: "
    print "(A, I0, A)", achar(9), size(s3, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s3%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s3%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s3%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "ROSENBROCK: "
    print "(A, I0, A)", achar(9), size(s4, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s4%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s4%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s4%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s4%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s4%statistics%matrix_factorizations, " Matrix Factorizations"

    ! Now, implement a PI controller and check its effect.  This will likely
    ! increase the number of steps (loss of efficiency), but if there were
    ! any stability issues, stability will likely improve.  Stability is likely
    ! not relevant on this problem, but it's here for illustration purposes.
    call integrator_2%clear_buffer()
    call integrator_2%set_step_size_control_parameter(0.1d0)
    call integrator_2%solve(mdl, t, ic)
    s2a = integrator_2%get_ode_results()
    print "(A)", "RUNGE_KUTTA_45 w/ PI Controller: "
    print "(A, I0, A)", achar(9), size(s2a, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s2a%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s2a%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s2a%statistics%function_evaluations, " Function Evaluations"

    ! Kennedy-Carpenter Integrators
    print "(A)", "KC4: "
    print "(A, I0, A)", achar(9), size(s5, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s5%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s5%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s5%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s5%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s5%statistics%matrix_factorizations, " Matrix Factorizations"

    print "(A)", "KC5: "
    print "(A, I0, A)", achar(9), size(s6, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s6%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s6%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s6%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s6%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s6%statistics%matrix_factorizations, " Matrix Factorizations"

    ! Additional Integrators
    print "(A)", "TSITOURAS 4/5: "
    print "(A, I0, A)", achar(9), size(s7, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s7%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s7%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s7%statistics%function_evaluations, " Function Evaluations"

    print "(A)", "BDF: "
    print "(A, I0, A)", achar(9), size(s8, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s8%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s8%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s8%statistics%jacobian_evaluations, " Jacobian Evaluations"
    print "(A, I0, A)", achar(9), s8%statistics%function_evaluations, " Function Evaluations"
    print "(A, I0, A)", achar(9), s8%statistics%matrix_factorizations, " Matrix Factorizations"

    print "(A)", "ADAMS: "
    print "(A, I0, A)", achar(9), size(s9, 1), " Solution Points"
    print "(A, I0, A)", achar(9), s9%statistics%successful_steps, " Successful Steps"
    print "(A, I0, A)", achar(9), s9%statistics%failed_steps, " Failed Steps"
    print "(A, I0, A)", achar(9), s9%statistics%function_evaluations, " Function Evaluations"
end program
```
```txt
RUNGE_KUTTA_23:
    2465 Solution Points
    2464 Successful Steps
    26 Failed Steps
    7471 Function Evaluations
RUNGE_KUTTA_45:
    583 Solution Points
    582 Successful Steps
    87 Failed Steps
    4015 Function Evaluations
RUNGE_KUTTA_853:
    925 Solution Points
    924 Successful Steps
    84 Failed Steps
    12097 Function Evaluations
ROSENBROCK:
    1185 Solution Points
    1184 Successful Steps
    94 Failed Steps
    1184 Jacobian Evaluations
    11221 Function Evaluations
    1278 Matrix Factorizations
RUNGE_KUTTA_45 w/ PI Controller:
    1107 Solution Points
    1106 Successful Steps
    0 Failed Steps
    6637 Function Evaluations
KC4:
    483 Solution Points
    482 Successful Steps
    145 Failed Steps
    6244 Jacobian Evaluations
    31728 Function Evaluations
    6244 Matrix Factorizations
KC5:
    245 Solution Points
    244 Successful Steps
    73 Failed Steps
    4840 Jacobian Evaluations
    24042 Function Evaluations
    4840 Matrix Factorizations
TSITOURAS 4/5:
    522 Solution Points
    521 Successful Steps
    102 Failed Steps
    4362 Function Evaluations
BDF:
    1584 Solution Points
    1583 Successful Steps
    59 Failed Steps
    3276 Jacobian Evaluations
    16331 Function Evaluations
    3276 Matrix Factorizations
ADAMS:
    1133 Solution Points
    1132 Successful Steps
    70 Failed Steps
    2336 Function Evaluations
```

Because the integration range is supplied as just its two end points, each solver stores the initial condition and one point per accepted step.  The solution-point count is therefore one larger than `statistics%successful_steps` in this example.  With dense output, the stored points instead follow the requested output grid, so their count need not reflect the number of internal integration steps at all.  The statistics report the solver's work independently of that grid: successful steps show how many steps advanced the integration, failed steps expose retries that produced no solution point, and Jacobian evaluations reveal additional work performed by implicit methods.

Function evaluations and matrix factorizations make the cost comparison more informative because the work performed per step differs considerably from one method to the next.  Explicit Runge-Kutta methods evaluate the model at multiple stages; implicit methods may also form Jacobians, factor matrices, and iterate to convergence.  Comparing `function_evaluations` alongside step counts reveals how many model evaluations a method required, rather than inferring that work from its accepted steps.  The separate `matrix_factorizations` counter exposes work required to prepare linear systems for solution; it need not equal the number of Jacobian evaluations or successful steps.  These operations have different costs depending on the model, system size, and method, and the statistics do not measure elapsed time.  Fewer evaluations or factorizations therefore do not necessarily imply a shorter wall-clock time.

For example, `rosenbrock` reports 1278 matrix factorizations but only 1184 Jacobian evaluations in this run.  The counters describe distinct operations rather than interchangeable estimates of implicit-solver work.  KC4, KC5, and BDF happen to report equal Jacobian-evaluation and factorization counts here; that equality should not be assumed for every method.

With $\mu = 5$ the Van der Pol oscillator is only mildly stiff, so the explicit integrators remain competitive.  The low-order `runge_kutta_23` requires by far the most steps because its step size is limited by accuracy rather than by stability; the higher-order `tsitouras_54` and `runge_kutta_45` traverse the same interval in roughly a fifth as many steps.  The Kennedy-Carpenter methods accept the fewest steps of any integrator here, and `kennedy_carpenter_5` accepts the fewest of all, but each of those steps carries the cost of the Newton iteration noted above.  The implicit methods may accept fewer steps while doing more work per step than the explicit methods.

The two multi-step methods make a useful counterpoint, because this problem is close to their worst case.  Van der Pol at $\mu = 5$ is a relaxation oscillator: long, smooth excursions punctuated by abrupt transitions.  A multi-step method earns its efficiency by carrying a history of accepted points and raising its order as that history proves reliable, and both `bdf` and `adams` restart at first order whenever the step is disturbed.  Repeated sharp transitions therefore deny them the sustained high order at which they excel, which is why `bdf` records the largest step count of any implicit method here.  Their advantage appears instead on smooth problems held to tight tolerances, where high order is sustainable; the order-selection study behind the `adams` entry above is a better illustration of that regime than this example is.

Step counts are especially misleading for `adams`, which uses neither a Jacobian nor a linear solve.  It takes 1132 successful steps and 2336 function evaluations here, whereas `kennedy_carpenter_5` takes only 244 successful steps but requires 24042 function evaluations, 4840 Jacobian evaluations, and 4840 matrix factorizations.  The function-evaluation count makes this difference visible without assuming a fixed number of evaluations per step.  Comparing failed steps as well as successful steps exposes retries, while Jacobian evaluations and matrix factorizations distinguish some of the extra work done by the implicit methods.  These counters complement one another, but a ranking by wall-clock time still requires timing the same problem at comparable accuracy.

The final result illustrates the PI step-size controller.  Applying it to `runge_kutta_45` raises the successful-step count from 582 to 1106, but reduces failed steps from 87 to zero for this problem.  Function evaluations nevertheless increase from 4015 to 6637.  The statistics make this tradeoff explicit: a smoother sequence of step sizes can avoid retries while still requiring more accepted steps and more model evaluations overall.  This can be valuable when the error estimate is noisy or when stability rather than accuracy is limiting the step, but fewer failed steps alone do not establish greater efficiency.  Stability is not a concern for this problem, so the controller is shown here purely for illustration.  PI control is disabled by default and must be requested explicitly.  Note that the controller applies to the single-step integrators; `rosenbrock` and the two multi-step methods use step-size estimators of their own.

### Differential-Algebraic Equations
This example illustrates the solution of a differential-algebraic equation (DAE) by means of a singular mass matrix.  The model is a simple pendulum expressed in Cartesian coordinates.  A mass $m$ swings from a pivot at the origin on a rigid, massless rod of length $L$, so the position of the mass $(x, y)$ must satisfy the constraint

$$
x^2 + y^2 = L^2.
$$

Introducing a Lagrange multiplier $\lambda$ to account for the constraint force, the equations of motion are

$$
\dot{x} = v_x, \qquad m \dot{v}_x = 2 \lambda x,
$$

$$
\dot{y} = v_y, \qquad m \dot{v}_y = 2 \lambda y - m g.
$$

Collecting the state into $z = \left[ x, v_x, y, v_y, \lambda \right]^T$, the system takes the mass matrix form $M \dot{z} = f(t, z)$ where

$$
M = \begin{bmatrix}
1 & 0 & 0 & 0 & 0 \\
0 & m & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & m & 0 \\
0 & 0 & 0 & 0 & 0
\end{bmatrix}.
$$

The zero in the final diagonal entry is what distinguishes this problem from an ordinary differential equation.  The multiplier has no time derivative of its own, so $M$ is singular and the last row states an algebraic relationship rather than a differential one.  That relationship comes from the constraint itself.  As written, the constraint involves only position, so it must be differentiated twice with respect to time before the multiplier appears in it.  Substituting the equations of motion into the result gives

$$
0 = 2 \lambda \left( x^2 + y^2 \right) + m \left( v_x^2 + v_y^2 - g y \right),
$$

which is the fifth equation of the system.  Because

$$
\frac{\partial}{\partial \lambda} \left[ 2 \lambda \left( x^2 + y^2 \right) + m \left( v_x^2 + v_y^2 - g y \right) \right] = 2 \left( x^2 + y^2 \right) = 2 L^2 \neq 0,
$$

this equation determines $\lambda$, and the problem is therefore of index 1.  The multiplier is not computed by the model.  It is carried as a state and the solver determines it, along with the four differential states, from the coupled system.  Two consequences are worth noting.  The initial conditions must be consistent, meaning they must satisfy both the constraint and the algebraic equation, and because the constraint has been differentiated the position constraint is only enforced to within the accuracy of the integration.  A slow drift in $x^2 + y^2 - L^2$ is expected over a long enough integration.

The model parameters are passed to the solver by means of the optional `args` argument.  A derived type is used to carry both quantities.
```fortran
type cartesian_pendulum_properties
    real(real64) :: mass
    real(real64) :: length
end type
```

The routine defining the system returns the four derivatives followed by the residual of the algebraic equation.  Notice that $\lambda$, carried as `x(5)`, is used rather than computed.
```fortran
subroutine cartesian_pendulum(t, x, dxdt, args)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt
    class(*), intent(inout), optional :: args

    ! Local Variables
    real(real64), parameter :: gc = 9.81d0
    real(real64) :: m, L

    ! Model Parameters
    select type (args)
    class is (cartesian_pendulum_properties)
        L = args%length
        m = args%mass
    end select

    ! The state vector is [x, dx/dt, y, dy/dt, lambda], where lambda is the
    ! Lagrange multiplier enforcing the constraint equation
    ! x**2 + y**2 = L**2.  The constraint force acts along the rod, and is
    ! given by lambda times the gradient of the constraint equation.
    !
    ! The constraint is differentiated twice with respect to time to reduce
    ! the system to index 1, which supplies the fifth equation below.  Notice,
    ! lambda is not computed here.  The fifth equation is written as a
    ! residual that the solver drives to zero, and the corresponding row of
    ! the mass matrix is zero.  It is that zero row which makes the system a
    ! DAE rather than an ODE.
    dxdt(1) = x(2)
    dxdt(2) = 2.0d0 * x(5) * x(1)
    dxdt(3) = x(4)
    dxdt(4) = 2.0d0 * x(5) * x(3) - m * gc
    dxdt(5) = 2.0d0 * x(5) * (x(1)**2 + x(3)**2) + &
        m * (x(2)**2 + x(4)**2 - gc * x(3))
end subroutine
```

The mass matrix routine returns $M$.  The matrix is constant for this model, so the solver can be told that it need only be evaluated once.
```fortran
subroutine cartesian_pendulum_mass_matrix(t, x, m, args)
    ! Arguments
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:,:) :: m
    class(*), intent(inout), optional :: args

    ! Parameters
    real(real64) :: mass, L

    ! Model Parameters
    select type (args)
    class is (cartesian_pendulum_properties)
        L = args%length
        mass = args%mass
    end select

    m = 0.0d0
    m(1,1) = 1.0d0
    m(2,2) = mass
    m(3,3) = 1.0d0
    m(4,4) = mass
    m(5,5) = 0.0d0  ! algebraic constraint equation
end subroutine
```

The mass matrix is supplied to the `ode_container` alongside the routine defining the system.  The pendulum is released from rest in a horizontal position at $x = L$, $y = 0$.  These conditions are consistent, as they satisfy the constraint, and with the mass at rest at $y = 0$ the algebraic equation reduces to $2 \lambda L^2 = 0$, so the multiplier begins at zero.
```fortran
program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Model Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: mass = 2.0d0

    ! Local Variables
    type(rosenbrock) :: integrator
    type(ode_container) :: mdl
    type(cartesian_pendulum_properties) :: args
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd

    ! Pass in the model parameters
    args%length = length
    args%mass = mass

    ! Define the model
    mdl%fcn => cartesian_pendulum
    mdl%mass_matrix => cartesian_pendulum_mass_matrix
    call mdl%set_is_mass_matrix_dependent(.false.)

    ! Compute the solution
    call integrator%solve( &
        mdl, &                                  ! model to solve
        [0.0d0, 1.0d1], &                       ! time bounds
        [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], & ! initial conditions - must satisfy algebraic constraints
        args = args &                           ! arguments to pass to the model
    )
    sol = integrator%get_solution()

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()
    call plt%set_use_y2_axis(.true.)
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call y2Axis%set_title("y(t)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pd1%define_data(sol(:,1), sol(:,2))
    call pd1%set_name("x(t)")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(sol(:,1), sol(:,4))
    call pd2%set_name("y(t)")
    call pd2%set_line_width(2.0)
    call pd2%set_draw_against_y2(.true.)
    call plt%push(pd2)

    call plt%draw()
end program
```
![](images/dae_results.png?raw=true)

The mass swings between $x = \pm L$ while $y$ remains at or below the pivot, which is the expected behavior for a pendulum released from rest in a horizontal position ($\theta = 0$).

The example uses `rosenbrock`, but any integrator that tolerates a singular mass matrix will serve.  Substituting `type(bdf)`, `type(kennedy_carpenter_4)`, or `type(kennedy_carpenter_5)` for the declared integrator type is the only change required, as every solver exposes the same `solve` and `get_solution` interface.

### Stabilized Index-1 GGL

Use `ggl_container` instead of `ode_container` for mechanical systems with
autonomous holonomic constraints $g(q)=0$. This is a **stabilized index-1
GGL-inspired variant**, not classical index-2 GGL. It replaces the original
position and velocity algebraic constraints with their stabilized evolution
equations. For $N$ coordinates and $M$ constraints, define

$$
z=[q_N,v_N,\lambda_M,\mu_M]^T,\qquad
G(q)=\frac{\partial g}{\partial q},\qquad
\mathcal M=diag(I_N,M_N,0_M,0_M).
$$

The container supplies the equations

$$
\begin{aligned}
\dot q &= v+G^T\mu,\\
M\dot v &= F(x,q,v)-G^T\lambda,\\
0 &= G(v+G^T\mu)+\alpha g(q),\\
0 &= DG(q)[v+G^T\mu]v+GM^{-1}(F-G^T\lambda)+\beta Gv.
\end{aligned}
$$

The force multiplier uses the sign convention $-G^T\lambda$. The last two
RHS blocks occupy the rows corresponding to the multiplier slots, but contain
no multiplier derivatives. The physical $M$ must be nonsingular; only the
augmented matrix is singular. Independent constraints and a symmetric
positive-definite physical mass guarantee the required regularity of $GG^T$
and $GM^{-1}G^T$. Singular physical or initialization systems are rejected.

These equations imply

$$
\dot g=-\alpha g,\qquad \frac{d(Gv)}{dx}=-\beta Gv.
$$

Both constraints remain zero in the exact solution if initially satisfied.
Positive gains damp initial or numerical violations; zero gains preserve
the corresponding initial residual rather than correcting it. Gains default
to one and have inverse-time units when $x$ is time. Large gains can introduce
additional stiffness. Stabilization is not projection: neither accepted steps
nor interpolated output are guaranteed to lie exactly on the constraint surface.
Explicitly time-dependent constraints are not supported.

The callback shapes differ from an ordinary `ode_container`:

| Callback | Input State | Output |
| --- | --- | --- |
| `fcn` (required) | `[q,v]`, length $2N$ | Applied force $F$, length $N$, excluding multiplier forces. |
| `mass_matrix` (optional) | $q$, length $N$ | Physical $N\times N$ mass; identity if absent. |
| `constraints` (required) | $q$, length $N$ | Position residual $g$, length $M$. |
| `constraint_jacobian` (required) | $q$, length $N$ | $G$, an $M\times N$ matrix. |
| `constraint_curvature` (required) | $q$, direction $\dot q$, and physical velocity $v$, each length $N$ | $DG(q)[\dot q]v$, length $M$; zeros for linear constraints. |
| `jacobian` (optional) | Full augmented state $z$ | Jacobian of the complete augmented RHS, of size $(2N+2M)\times(2N+2M)$. |

All callbacks receive the independent variable and optional user `args`.
The constraint callbacks must describe autonomous geometry, accept off-manifold
inputs, and fill every output element. There is no finite-difference fallback
for $G$ or curvature. If the full RHS Jacobian is omitted, the inherited
finite-difference implementation evaluates the augmented RHS automatically.
Set `set_is_mass_matrix_dependent(.true.)` when the physical mass varies with
time or coordinates; otherwise solvers can cache the augmented mass per solve.

`initialize(N,M,alpha,beta)` defines dimensions and gains; `get_state_count()`
returns $2N+2M$. `validate_initial_conditions(x,z,tolerance,args)` checks both
position and velocity residuals against an absolute tolerance without changing
the state. Choose constraint scaling that makes this tolerance meaningful for
both residuals. `initialize_multipliers(x,z,args)` solves for $\mu$ and $\lambda$
without changing $q$ or $v$. It makes the stabilized algebraic rows consistent,
but does not repair physically inconsistent position or velocity inputs.

Here is a complete pendulum example with $g(q)=q^Tq-L^2$, $G=2q^T$, and
$DG(q)[\dot q]v=2\dot q^Tv$. Its state ordering differs from the preceding
acceleration-level pendulum example: `[qx,qy,vx,vy,lambda,mu]`.
The complete source is available in
[stabilized_pendulum_example.f90](examples/stabilized_pendulum_example.f90).
It uses FPLOT to plot the Cartesian position and pendulum angle over 10 seconds.

```fortran
program stabilized_pendulum
    use iso_fortran_env, only : real64
    use diffeq, only : ggl_container, bdf
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Parameters
    real(real64), parameter :: length = 1.5d0
    real(real64), parameter :: body_mass = 2.0d0

    ! Local Variables
    type(ggl_container) :: model
    type(bdf) :: solver
    real(real64), dimension(6) :: initial_state
    real(real64), allocatable, dimension(:) :: theta
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: pltxy, plttheta
    class(legend), pointer :: lgnd

    ! Set up the model and constraints
    call model%initialize(2, 1, alpha = 2.0d0, beta = 3.0d0)
    model%fcn => force
    model%mass_matrix => physical_mass
    model%constraints => position_constraint
    model%constraint_jacobian => constraint_gradient
    model%constraint_curvature => constraint_curvature
    call model%set_is_mass_matrix_dependent(.false.)

    ! Define the initial conditions and validate against the constraints.  The
    ! initial conditions must satisfy the constraints.
    initial_state = [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    if (.not.model%validate_initial_conditions(0.0d0, initial_state, 1.0d-10)) &
        error stop "Inconsistent physical initial conditions"
    call model%initialize_multipliers(0.0d0, initial_state)

    ! Set up the solver and adjust tolerances.  Utilize the same tolerance for
    ! every equation.
    call solver%set_absolute_tolerance(1.0d-9)
    call solver%set_relative_tolerance(1.0d-9)

    ! Solve the model and retrieve the results
    call solver%solve(model, [0.0d0, 1.0d1], initial_state)
    sol = solver%get_solution()

    ! Compute the angular representation of the solution
    theta = atan2(sol(:,3), sol(:,2)) * 1.8d2 / pi

! ------------------------------------------------------------------------------
    ! Set up the plot
    call plt%initialize(1, 2, width = 1200)

    ! Plot the X-Y results
    call pltxy%initialize()
    lgnd => pltxy%get_legend()
    call pltxy%set_x_axis_title("t")
    call pltxy%set_y_axis_title("x(t), y(t)")
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)

    call pltxy%push(sol(:,1), sol(:,2), name = "x(t)")
    call pltxy%push(sol(:,1), sol(:,3), name = "y(t)")

    call plt%set(1, 1, pltxy)

    ! Plot the angular results
    call plttheta%initialize()
    call plttheta%set_x_axis_title("t")
    call plttheta%set_y_axis_title("{/Symbol q} [deg]")
    call plttheta%push(sol(:,1), theta)
    
    call plt%set(1, 2, plttheta)

    call plt%draw()
contains
    subroutine force(x, y, applied_force, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: applied_force
        class(*), intent(inout), optional :: args
        applied_force = [0.0d0, -body_mass * 9.81d0]
    end subroutine

    subroutine physical_mass(x, q, mass, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: mass
        class(*), intent(inout), optional :: args
        mass = 0.0d0
        mass(1,1) = body_mass
        mass(2,2) = body_mass
    end subroutine

    subroutine position_constraint(x, q, g, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(inout), optional :: args
        g(1) = dot_product(q, q) - length**2
    end subroutine

    subroutine constraint_gradient(x, q, gradient, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: gradient
        class(*), intent(inout), optional :: args
        gradient(1,:) = 2.0d0 * q
    end subroutine

    subroutine constraint_curvature(x, q, direction, v, curvature, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q, direction, v
        real(real64), intent(out), dimension(:) :: curvature
        class(*), intent(inout), optional :: args
        curvature(1) = 2.0d0 * dot_product(direction, v)
    end subroutine
end program
```

Use BDF, Rosenbrock, or either Kennedy-Carpenter solver with this container.
Adams and explicit Runge-Kutta methods are unsuitable for its singular augmented
matrix. The regression suite checks constraint preservation and exponential
decay for all four supported implicit solvers, including dense output and a
time-dependent physical mass.

The stabilized formulation produces the following response. The left plot
shows the horizontal coordinate oscillating between $-L$ and $L$, while the
vertical coordinate stays between $-L$ and zero. The mass passes its lowest
point twice per complete horizontal oscillation, so the vertical motion repeats
twice as often. The right plot shows
$\theta$, measured from the positive horizontal
axis: the release point is $0^\circ$, the lowest point is $-90^\circ$, and
the opposite horizontal turning point is $-180^\circ$.

![Stabilized GGL pendulum Cartesian coordinates and angle over 10 seconds](images/stabilized_pendulum_example.png?raw=true)

For comparison, the preceding acceleration-level formulation produces:

![Acceleration-level pendulum Cartesian coordinates and angle over 10 seconds](images/dae_results.png?raw=true)

At the resolution of these plots, the two formulations show closely matching
amplitudes and oscillation timing, with no obvious decay of the swing amplitude
over the displayed interval. This is the expected agreement for the same
undamped pendulum released from consistent initial conditions. The GGL gains
are not physical damping coefficients: they act on constraint violations,
not on the intended pendulum motion.

The distinction is how departures from the constraint surface are treated.
The earlier model enforces the acceleration-level constraint, but does not
explicitly damp position or velocity constraint errors. The stabilized variant
instead gives the residuals
$r_p=q_x^2+q_y^2-L^2$ and $r_v=2(q_xv_x+q_yv_y)$ the evolution laws
$\dot r_p=-2r_p$ and $\dot r_v=-3r_v$ in the exact solution for the gains
used here. These relationships provide restoring behavior for violations, without
projecting the numerical state onto the constraint surface.

The motion plots alone cannot quantify that benefit or establish that one
formulation is more accurate or efficient. A numerical comparison should also
plot $r_p$ and $r_v$, measure energy drift, and compare solver work at matched
tolerances. Small constraint violations or accumulated phase errors can remain
invisible on the position and angle scales shown here.

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
9. Schoeder, Svenja & Ulbrich, Heinz & Schindler, Thorsten. (2013). Discussion of the Gear–Gupta–Leimkuhler method for impacting mechanical systems. Multibody System Dynamics. 31. 10.1007/s11044-013-9370-y. 
10. Pöll, C. & Hafner, Irene. (2015). Index Reduction and Regularisation Methods for Multibody Systems. IFAC-PapersOnLine. 48. 306-311. 10.1016/j.ifacol.2015.05.150. 
11. Tsitouras, Charalampos. (2011). Runge–Kutta pairs of order 5(4) satisfying only the first column simplifying assumption. Computers & Mathematics with Applications. 62. 770-775. 10.1016/j.camwa.2011.06.002. 
12. M. H. Carpenter, C. A. Kennedy, and J. M. Derlaga, “Intrastep, Stage-Value Predictors for Diagonally-Implicit Runge–Kutta Methods,” NASA Langley Research Center, Presentation, 2023.
13. Gear, C. W. (1971). Numerical initial value problems in ordinary differential equations. Prentice-Hall.
14. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving ordinary differential equations I: Nonstiff problems (2nd revised ed.). Springer.
15. Hairer, E., & Wanner, G. (1996). Solving ordinary differential equations II: Stiff and differential-algebraic problems (2nd revised ed.). Springer.
16. Gear, C. W., Leimkuhler, B., & Gupta, G. K. (1985). Automatic integration of Euler-Lagrange equations with constraints. Journal of Computational and Applied Mathematics, 12-13, 77-90. https://doi.org/10.1016/0377-0427(85)90008-1.
17. Baumgarte, J. (1972). Stabilization of constraints and integrals of motion in dynamical systems. Computer Methods in Applied Mechanics and Engineering, 1(1), 1-16. https://doi.org/10.1016/0045-7825(72)90018-7.
