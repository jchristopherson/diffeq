## Multistep Integrators

Multistep integrators advance the solution using a history of previously accepted points. The abstract `multistep_integrator` type in [diffeq_multistep.f90](diffeq_multistep.f90) provides history management, adaptive step-size and order selection, dense output, and the common `solve` and `get_solution` interface. Both methods start at order one and increase the order as sufficient history becomes available.

| Integrator | Method | Order | Intended problems | Source |
|---|---|---|---|---|
| `adams` | Adams--Bashforth--Moulton PECE | 1--12 | Non-stiff ODEs; nonsingular mass matrices | [diffeq_pece.f90](diffeq_pece.f90) |
| `bdf` | Backward differentiation formula | 1--5 | Stiff ODEs and index-1 mass-matrix DAEs | [diffeq_bdf.f90](diffeq_bdf.f90) |

### Adams Predictor--Corrector

`adams` uses an Adams--Bashforth prediction, evaluates the derivative at the predicted state, applies one Adams--Moulton correction, and evaluates the corrected derivative. The predictor--corrector difference estimates the local error. This PECE cycle does not require a Jacobian or a nonlinear solve; its single correction makes it best suited to non-stiff problems. Stiff problems may force small steps for stability.

The method supports a nonsingular mass matrix, which is factored to obtain the derivative. A singular mass matrix requires a DAE-capable method such as `bdf`.

### Backward Differentiation Formula

`bdf` approximates the derivative using a polynomial through the new point and previous solution points. Its coefficients use the actual history spacing, so the step size can vary. The implicit corrector is solved by Newton iteration. For $M(x,y)y' = f(x,y)$, the iteration matrix is approximately

$$
A = \frac{\alpha_0}{h}M - J, \qquad J = \frac{\partial f}{\partial y}.
$$

When no mass matrix is supplied, $M$ is the identity. The iteration matrix does not require $M$ to be inverted, so `bdf` can solve regular index-1 DAEs with singular mass matrices. Initial conditions should satisfy the algebraic constraints. For a state-dependent mass matrix, the iteration matrix neglects derivatives of the mass matrix; this can affect Newton convergence.

The dense path can use a supplied Jacobian, automatic differentiation, or finite differences. Supplying a sparse Jacobian selects the sparse Newton path, which uses CSR matrices and an ILUT-preconditioned GMRES solve. Sparse Jacobians must currently be provided by the user; automatic sparse Jacobian construction is not yet available. A mass matrix used with the sparse path must also be supplied as a sparse matrix.

### Choosing and Configuring an Integrator

- **Non-stiff systems:** use `adams` when a variable-order predictor--corrector is appropriate.
- **Stiff systems or index-1 DAEs:** use `bdf`; provide a sparse Jacobian for large sparse systems.

Both integrators expose `set_absolute_tolerance`, `set_relative_tolerance`, and `set_maximum_order`. With two output times, `solve` stores accepted integration steps; with more than two, it interpolates to the requested times. The first solution column contains the independent variable, followed by the state variables. Call `clear_buffer` before a separate solve when previous output should not be retained.

For a broader comparison, see [Choosing an Integrator](../../README.md#choosing-an-integrator) in the main README.

### Examples

The examples below are located in the top-level [examples/](../../examples) folder and use fplot to display their results.

#### Sparse BDF: Heat Equation

[sparse_bdf_example.f90](../../examples/sparse_bdf_example.f90) discretizes the one-dimensional heat equation with zero boundary temperatures. Its 1,000 interior points yield a stiff system with a tridiagonal Jacobian. The example supplies that Jacobian as a CSR matrix, solves at five output times, and compares the numerical result with the exact decay of the discrete system's eigenmodes.

```fortran
type(bdf) :: integrator
type(ode_container) :: model
real(real64), allocatable :: solution(:,:)

model%fcn => heat_rhs
model%sparse_jacobian => heat_jacobian
call integrator%solve(model, times, initial)
solution = integrator%get_solution()
```

The callbacks and the `times` and `initial` arrays are defined in the example. The plotted curves show the temperature profile at each requested time.

![Sparse BDF heat equation profiles](../../images/sparse_bdf_example.png?raw=true)

#### BDF: Automatic, Finite-Difference, and Analytical Jacobians

[ad_example.f90](../../examples/ad_example.f90) solves the Van der Pol equation with $\mu = 5$ using `bdf`. It compares a dual-number right-hand side for automatic differentiation (AD), a real-valued right-hand side with a finite-difference Jacobian (FD), and an explicitly supplied analytical Jacobian. The example also reports timing and solver statistics.

```fortran
type(ode_container) :: modelAD, modelFD, model
type(bdf) :: integrator
type(ode_results) :: rstAD, rstFD, rst

modelAD%fcn_ad => vanderpol_ad
call integrator%solve(modelAD, x, ic)
rstAD = integrator%get_ode_results()

modelFD%fcn => vanderpol
call integrator%clear_buffer()
call integrator%solve(modelFD, x, ic)
rstFD = integrator%get_ode_results()

model%fcn => vanderpol
model%jacobian => vanderpol_analytical_jacobian
call integrator%clear_buffer()
call integrator%solve(model, x, ic)
rst = integrator%get_ode_results()
```

The figure compares the solutions and their differences from the analytical-Jacobian run; the latter is a numerical reference, not an exact solution of the ODE.

![BDF solutions and differences for three Jacobian approaches](../../images/ad_example.png?raw=true)

#### Adams PECE: Damped Oscillator

[adams_example.f90](../../examples/adams_example.f90) solves a non-stiff damped oscillator with initial displacement one and velocity zero. It compares displacement and velocity with the exact solution, and reports the accepted-step count and highest order used.

```fortran
type(adams) :: integrator
type(ode_container) :: model
real(real64), allocatable :: solution(:,:)

model%fcn => oscillator_rhs
call integrator%set_absolute_tolerance(1.0d-8)
call integrator%set_relative_tolerance(1.0d-8)
call integrator%solve(model, times, [1.0d0, 0.0d0])
solution = integrator%get_solution()
```

![Adams solution and exact damped oscillator response](../../images/adams_example.png?raw=true)