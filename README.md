# diffeq
A modern Fortran library providing an object-oriented approach to solving and exploring ordinary differential equations.

## Build Status
[![CMake](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/diffeq/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/diffeq/workflows/fpm/badge.svg)](https://github.com/jchristopherson/diffeq/actions)

## Documentation
- [Fortran API Documentation](https://jchristopherson.github.io/diffeq/).
- [C API Documentation](https://jchristopherson.github.io/diffeq/c-api/).

## Capabilities
| Capability | Summary |
| --- | --- |
| Ordinary differential equations (ODEs) | Solve initial-value problems with adaptive step-size control and configurable absolute and relative tolerances. |
| Mass-matrix systems | Formulate problems as $M y' = f(x,y)$, including systems with state-dependent mass matrices. |
| Index-1 differential-algebraic equations (DAEs) | Solve supported singular-mass-matrix systems with consistent initial conditions. |
| Holonomically constrained mechanics | Use a stabilized index-1 GGL-inspired formulation with position and velocity constraint residual damping. |
| Dense output | Request a specified output grid or evaluate the last prepared accepted step through a common solver interface. |
| Event detection | Stop integration at zero crossings of user-defined event functions, then restart with updated initial conditions if needed. |
| Solution diagnostics | Retrieve step counts, function and Jacobian evaluations, and matrix-factorization statistics. |
| C interoperability | Integrate from C through the library's C API. |

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
and a consistent-derivative calculation. Their stage solves use a simplified
Newton iteration that reuses a single factorization of $M - \gamma h J$ for
every stage, keeps the Jacobian across steps while the iteration converges
quickly, and refreshes it only when convergence slows or fails. `bdf` factors
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
| `model%fcn` | Defines $f(x,y)$ with real-valued arguments. Assign either this callback or `model%fcn_ad`, not both. |
| `model%fcn_ad` | Defines the RHS using autodiff `dual` arguments, enabling automatic differentiation of its state and time dependence. |
| `model%set_thread_count(n)` | Sets the maximum OpenMP thread count for AD or finite-difference Jacobian calculations (default 1). The RHS callback must be thread-safe when using multiple threads. |
| `model%jacobian` | Supplies $J = \partial f / \partial y$. Takes precedence over AD and finite differences; otherwise the selected RHS callback determines the differentiation method. |
| `model%mass_matrix` | Supplies $M(x,y)$ for a system written as $M y' = f(x,y)$. It is supported by the Rosenbrock, Kennedy-Carpenter, BDF, and Adams solvers. Only Adams requires $M$ to be nonsingular. |
| `model%event_function` | Supplies one or more event values; DIFFEQ locates zero crossings and stops the solve at a detected event. Set the number of event values with `model%set_event_equation_count(count)`. |
| `integrator%solve(model, x, y0, store_iv)` | Integrates from `x(1)` to `x(size(x))` using initial state `y0`. The optional `store_iv` defaults to true; set it to false when restarting at an event already present in the result buffer. |
| `integrator%get_solution()` | Returns an $N \times (n+1)$ array with the independent variable in column 1 and state values in columns 2 through $n+1$. |
| `integrator%get_ode_results()` | Returns an `ode_results` container that provides the solution and statistics on integrator performance, including step counts and function, Jacobian, and matrix-factorization evaluations. |
| `integrator%set_dense_output_enabled(enabled)` | Opts into dense-output preparation on every accepted step. The default is false. |
| `integrator%has_dense_output()` | Tests whether a prepared accepted-step interpolant is available. |
| `integrator%get_dense_interval(x_start, x_end)` | Returns the prepared step's endpoints in integration order, including for backward integration. |
| `integrator%evaluate_dense(x, y)` | Evaluates the prepared step at `x` into the state-sized array `y`, without advancing integration or evaluating the ODE. |
| `set_absolute_tolerance` / `set_relative_tolerance` | Control the local error scale, approximately $\mathrm{atol} + \mathrm{rtol}|y|$. |
| `set_step_size_control_parameter` | Sets the PI controller parameter used by the Kennedy-Carpenter solvers. Rosenbrock and the multi-step solvers use their own step-size estimators. |
| `set_maximum_order` / `get_maximum_order` | Bound the order of a multi-step solver. A request above what the method supports is clamped rather than rejected. |
| `set_newton_step_limit` / `set_newton_tolerance` | Govern the Newton iteration used by `bdf` and the Kennedy-Carpenter solvers. For the Kennedy-Carpenter solvers the defaults are 12 iterations and a tolerance of $10^{-8}$; the effective tolerance is further capped at 1% of the smallest integration tolerance so that the Newton error stays below the local error estimate. |
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

### Automatic Differentiation

Fortran callers can implement their RHS once with `dual` arguments from
[autodiff](https://github.com/jchristopherson/autodiff) and assign
it to `model%fcn_ad` instead of `model%fcn`. The integrator still uses real-valued
states, solutions, tolerances, mass matrices, and event callbacks:

```fortran
program ad_example
    use diffeq
    use autodiff
    implicit none
    type(ode_container) :: model
    type(rosenbrock) :: solver
    real(real64), allocatable :: solution(:,:)

    model%fcn_ad => rhs
    call solver%solve(model, [0.0d0, 1.0d0], [1.0d0])
    solution = solver%get_solution()
    print *, solution(size(solution,1),:)
contains
    subroutine rhs(x, y, dydx, args)
        !! Both the independent variable and state use AD arithmetic.
        type(dual), intent(in) :: x, y(:)
        type(dual), intent(out) :: dydx(:)
        class(*), intent(in), optional :: args
        dydx(1) = -2.0d0 * y(1) + exp(-x)
    end subroutine
end program
```

The callback must assign every output element using AD expressions or constant
assignment (which clears the derivative). Assigning only an output's `%x` is
insufficient. The type stores a value `%x` and one directional derivative `%dx`;
scalar and array expressions support arithmetic, comparisons, and common real 
mathematical intrinsics. Comparisons use values only. At nonsmooth points, `abs`
uses derivative zero at zero, binary `min`/`max` select the first argument on 
ties, and `sign` treats the sign selection as constant. Underlying 
real-intrinsic domain restrictions still apply; AD does not make 
nondifferentiable expressions differentiable.

Exactly one RHS representation is required. Assigning both callbacks raises
`DIFFEQ_INVALID_OPERATION_ERROR`, even if an explicit Jacobian is supplied.
An explicit `model%jacobian` always takes precedence; otherwise `fcn_ad` selects
AD and `fcn` retains the existing forward-difference behavior. Existing real
callbacks and the C API need no changes.

Ordinary RHS evaluation uses zero derivative seeds and extracts values. For
an N-state Jacobian, the container holds time fixed and seeds one state per
column, requiring N callback evaluations. Rosenbrock and Kennedy-Carpenter DAE
initialization also use a time seed to obtain partial $\partial f/\partial x$
in one callback evaluation. The container uses autodiff's reusable workspace
and `jacobian` routine.

OpenMP is required by both the CMake and fpm builds. AD and finite-difference
Jacobians are serial by default. Call `model%set_thread_count(n)` with a
positive integer to evaluate columns in parallel; `model%get_thread_count()`
returns the setting. The same setting applies to whichever differentiation
method the model uses. At most N threads are requested for an N-state
Jacobian, and OpenMP runtime limits may reduce the actual team size. With a
thread count of one, a serial code path that makes no OpenMP calls is used.
For finite differences, the unperturbed RHS is evaluated once, after which
each thread perturbs its own copy of the state to compute its columns.
Ordinary RHS and time-derivative evaluations remain serial. An explicit
Jacobian is unaffected by this setting. Parallelism is usually useful for
larger, more expensive models; thread overhead can outweigh its benefit for
small systems.

Optional `args` are `intent(in)` and forwarded unchanged. Repeated evaluations
at the same numerical inputs must describe the same function regardless of
derivative seeds. When using multiple threads, the RHS callback (`fcn` or
`fcn_ad`) must be thread-safe, including any code it calls: avoid
unsynchronized logging, counters, mutable shared data, or other side effects.
Do not concurrently use the same container for separate solves or Jacobian
requests.

`statistics%function_evaluations` counts actual RHS callback invocations,
including N calls per AD Jacobian, N+1 per finite-difference Jacobian, and one
per partial time derivative. An explicit Jacobian contributes one
`jacobian_evaluations` entry and no library-managed RHS calls. Calls made
internally by user-written Jacobian code are not tracked. These are counts of
evaluations, not estimates of equal computational cost.

For direct container use, `compute_jacobian` and `compute_rhs` accept real
inputs for either representation. `get_jacobian_function_evaluations()` reports
the RHS calls made by the last `compute_jacobian`, initially zero.
`compute_time_derivative(x, y, f, dfdx, args)` computes partial time dependence;
for a real callback, supply the already evaluated `f(x,y)` for forward
differencing. Container extensions overriding Jacobian evaluation should
override the evaluation-count query when they perform additional RHS work.

This feature differentiates the RHS, not the entire mass-matrix residual.
State/time derivatives of mass matrices retain their existing finite-difference
treatment. `ggl_container` does not support `fcn_ad`: its specialized force
callback and real-valued constraint assembly require a separate AD design.
It rejects `fcn_ad` explicitly; continue using its existing `fcn` and optional
full augmented-system Jacobian.

### Local Dense Output (Fortran)

All integrators expose the same interface for evaluating the last prepared
accepted step, including through a `class(ode_integrator)` reference:

```fortran
call solver%set_dense_output_enabled(.true.)
call solver%solve(model, [0.0d0, 1.0d0], initial_state, args)
if (solver%has_dense_output()) then
    call solver%get_dense_interval(x_start, x_end)
    call solver%evaluate_dense(x_start + 0.5d0 * (x_end - x_start), y)
end if
```

Declare `y` with one element per state variable. The interval includes both
endpoints, and its endpoints follow integration order: `x_start > x_end` for
backward integration. This is a local interpolant, not a continuous solution
archive for the whole integration.

Explicit preparation is disabled by default to preserve the cost of existing
two-point solves. An output grid with more than two points, or interpolation
back from an overshoot, prepares dense output automatically. The
`get_dense_output_enabled` setting reports the explicit opt-in, whereas
`has_dense_output` reports actual availability. In particular, event detection
must not infer availability solely from the number of requested output points.
RK853 preparation requires three additional RHS evaluations per prepared step;
subsequent `evaluate_dense` calls require no RHS evaluations for any solver.

The interpolant is invalidated when another step is attempted, a new solve
starts, `clear_buffer` is called, the explicit preparation setting is set, or
multistep history is reset. It remains available after `solve` returns if the
final accepted step was prepared. Rejected steps do not publish dense output.
Multistep dense output caches the accepted step's polynomial, independently of
subsequent order selection. Interpolant accuracy remains method-specific.

Calling `evaluate_dense` or `get_dense_interval` without a prepared step reports
`DIFFEQ_INVALID_OPERATION_ERROR`. A wrongly sized output array reports
`DIFFEQ_ARRAY_SIZE_ERROR`; nonfinite coordinates or coordinates outside the
prepared interval report `DIFFEQ_INVALID_INPUT_ERROR`. No extrapolation is
performed by this interface.

Existing `interpolate` bindings retain their signatures and behavior. Custom
types derived directly from `ode_integrator` must implement the new
`compute_dense_output(this, x, y)` hook with `this` declared `intent(in)`.
Solver implementations publish a prepared step with
`set_dense_interval(x_start, x_end, state_count)` and invalidate it before
altering interpolation data with `invalidate_dense_output`. These are
implementation hooks; application code should use `evaluate_dense`.

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

The C API tests can be run separately with
`ctest --test-dir build --output-on-failure -R '^diffeq_c_'`. They are organized as:

- [diffeq_c_api_tests.c](test/diffeq_c_api_tests.c): default settings, result
  access, matrix layout, error codes, and file output.
- [diffeq_c_integrator_tests.c](test/diffeq_c_integrator_tests.c): table-driven
  accuracy, forward/reverse integration, overshoot-setting, callback statistics,
  and invalid-input tests shared by all C integrators. Each method has its own
  CTest entry; add new methods to both the function table and the CMake test list.
- [diffeq_c_rosenbrock_tests.c](test/diffeq_c_rosenbrock_tests.c): all combinations
  of optional analytical Jacobian and mass-matrix callbacks, including
  column-major matrix indexing and solver-specific statistics.

The example programs are enabled separately:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=ON -DBUILD_DIFFEQ_EXAMPLES=ON
cmake --build build --config Release
```

Install the library, generated Fortran module files, the C API header
(`diffeq.h`), CMake package files, and pkg-config metadata to a local prefix
with:

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

Each family of integrators has its own README describing the available solvers and walking through the related examples. The example source files are located in the [examples/](examples) folder, and the resulting plots are located in the [images/](images) folder.

| Topic | Description |
|-------|-------------|
| [Single-Step Integrators](src/singlestep/README.md) | Explicit Runge-Kutta, Rosenbrock, and Kennedy-Carpenter ESDIRK integrators, including mass-matrix DAEs and event detection. |
| [Multistep Integrators](src/multistep/README.md) | The variable-order BDF integrator, including its sparse Jacobian path, and the Adams PECE integrator. |
| [Index Reduction (DAEs)](src/dae/README.md) | The stabilized index-1 GGL formulation via `ggl_container` for constrained mechanical systems, with both dense and sparse examples. |
| [C API](src/C/README.md) | The C interface to the integrators, including callbacks, settings, results, and C examples. |

## Comparing Integrators
Here's another example comparing the behavior of several integrators on the same Van der Pol problem illustrated in the previous example.  Each integrator is used in the same manner; only its declared type changes.  Calling `get_ode_results()` retrieves an `ode_results` container holding both the solution and an `ode_results_statistics` member named `statistics`.  This provides a more informative comparison than the number of stored solution points alone:

- `statistics%successful_steps`: the number of successful integration steps.
- `statistics%failed_steps`: the number of unsuccessful step attempts.
- `statistics%jacobian_evaluations`: the number of Jacobian matrix evaluations, relevant to methods that use a Jacobian.
- `statistics%function_evaluations`: the number of RHS callback evaluations, including finite-difference and AD evaluations used to form Jacobians or partial time derivatives.
- `statistics%matrix_factorizations`: the number of matrix factorizations, relevant to methods that solve linear systems.

The following code solves the model with every integrator, retrieves each result, and repeats `runge_kutta_45` with a PI-type step-size controller. The print statements are omitted here; the complete source, including reporting, is in [comparison_example.f90](examples/comparison_example.f90).
```fortran
program example
    use iso_fortran_env, only : real64
    use diffeq
    use diffeq_models
    implicit none

    real(real64), parameter :: t(2) = [0.0d0, 5.0d1]
    real(real64), parameter :: ic(2) = [2.0d0, 0.0d0]
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

    mdl%fcn => vanderpol

    call integrator_1%solve(mdl, t, ic)
    s1 = integrator_1%get_ode_results()

    call integrator_2%solve(mdl, t, ic)
    s2 = integrator_2%get_ode_results()

    call integrator_3%solve(mdl, t, ic)
    s3 = integrator_3%get_ode_results()

    call integrator_4%solve(mdl, t, ic)
    s4 = integrator_4%get_ode_results()

    call integrator_5%solve(mdl, t, ic)
    s5 = integrator_5%get_ode_results()

    call integrator_6%solve(mdl, t, ic)
    s6 = integrator_6%get_ode_results()

    call integrator_7%solve(mdl, t, ic)
    s7 = integrator_7%get_ode_results()

    call integrator_8%solve(mdl, t, ic)
    s8 = integrator_8%get_ode_results()

    call integrator_9%solve(mdl, t, ic)
    s9 = integrator_9%get_ode_results()

    ! Repeat Runge-Kutta 4/5 with PI step-size control.
    call integrator_2%clear_buffer()
    call integrator_2%set_step_size_control_parameter(0.1d0)
    call integrator_2%solve(mdl, t, ic)
    s2a = integrator_2%get_ode_results()
end program
```

The complete example reports the following results:

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
    1101 Solution Points
    1100 Successful Steps
    37 Failed Steps
    1100 Jacobian Evaluations
    11260 Function Evaluations
    1137 Matrix Factorizations
RUNGE_KUTTA_45 w/ PI Controller:
    1107 Solution Points
    1106 Successful Steps
    0 Failed Steps
    6637 Function Evaluations
KC4:
    483 Solution Points
    482 Successful Steps
    145 Failed Steps
    207 Jacobian Evaluations
    11263 Function Evaluations
    559 Matrix Factorizations
KC5:
    246 Solution Points
    245 Successful Steps
    73 Failed Steps
    173 Jacobian Evaluations
    9278 Function Evaluations
    315 Matrix Factorizations
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

For example, `rosenbrock` reports 1137 matrix factorizations but only 1100 Jacobian evaluations in this run.  The counters describe distinct operations rather than interchangeable estimates of implicit-solver work.  BDF happens to report equal Jacobian-evaluation and factorization counts here, whereas the Kennedy-Carpenter methods reuse a single Jacobian across many steps and refactor only when the step size changes; that equality should not be assumed for every method.

With $\mu = 5$ the Van der Pol oscillator is only mildly stiff, so the explicit integrators remain competitive.  The low-order `runge_kutta_23` requires by far the most steps because its step size is limited by accuracy rather than by stability; the higher-order `tsitouras_54` and `runge_kutta_45` traverse the same interval in roughly a fifth as many steps.  The Kennedy-Carpenter methods accept the fewest steps of any integrator here, and `kennedy_carpenter_5` accepts the fewest of all, but each of those steps carries the cost of the Newton iteration noted above.  The implicit methods may accept fewer steps while doing more work per step than the explicit methods.

The two multi-step methods make a useful counterpoint, because this problem is close to their worst case.  Van der Pol at $\mu = 5$ is a relaxation oscillator: long, smooth excursions punctuated by abrupt transitions.  A multi-step method earns its efficiency by carrying a history of accepted points and raising its order as that history proves reliable, and both `bdf` and `adams` restart at first order whenever the step is disturbed.  Repeated sharp transitions therefore deny them the sustained high order at which they excel, which is why `bdf` records the largest step count of any implicit method here.  Their advantage appears instead on smooth problems held to tight tolerances, where high order is sustainable; the order-selection study behind the `adams` entry above is a better illustration of that regime than this example is.

Step counts are especially misleading for `adams`, which uses neither a Jacobian nor a linear solve.  It takes 1132 successful steps and 2336 function evaluations here, whereas `kennedy_carpenter_5` takes only 245 successful steps but requires 9278 function evaluations, 173 Jacobian evaluations, and 315 matrix factorizations.  The function-evaluation count makes this difference visible without assuming a fixed number of evaluations per step.  Comparing failed steps as well as successful steps exposes retries, while Jacobian evaluations and matrix factorizations distinguish some of the extra work done by the implicit methods.  These counters complement one another, but a ranking by wall-clock time still requires timing the same problem at comparable accuracy.

The final result illustrates the PI step-size controller.  Applying it to `runge_kutta_45` raises the successful-step count from 582 to 1106, but reduces failed steps from 87 to zero for this problem.  Function evaluations nevertheless increase from 4015 to 6637.  The statistics make this tradeoff explicit: a smoother sequence of step sizes can avoid retries while still requiring more accepted steps and more model evaluations overall.  This can be valuable when the error estimate is noisy or when stability rather than accuracy is limiting the step, but fewer failed steps alone do not establish greater efficiency.  Stability is not a concern for this problem, so the controller is shown here purely for illustration.  PI control is disabled by default and must be requested explicitly.  Note that the controller applies to the single-step integrators; `rosenbrock` and the two multi-step methods use step-size estimators of their own.

## External Libraries
Here is a list of external code libraries utilized by this library.  The CMake build script will include these dependencies automatically; however, it is highly recommended that an optimized BLAS and LAPACK already reside on your system for best performance (used by LINALG for linear algebra calculations).
- [linalg](https://github.com/jchristopherson/linalg)
- [autodiff](https://github.com/jchristopherson/autodiff) (dual numbers and parallel Jacobian calculations)
- OpenMP (Fortran compiler/runtime support, required)

CMake finds installed LINALG and autodiff packages or fetches reference
versions. An installed autodiff package must be built with
`AUTODIFF_ENABLE_OPENMP=ON` to support parallel Jacobians. The fetched version
has OpenMP enabled. OpenMP is a public link dependency and is discovered by
the installed DIFFEQ CMake package as well. The fpm manifest declares both
autodiff and the `openmp` metapackage, enabling the compiler's OpenMP flags.

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
