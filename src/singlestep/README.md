## Single-Step Integrators

Single-step integrators advance the solution from one point to the next using only the current state. Each integrator derives from the abstract `single_step_integrator` type defined in [diffeq_singlestep.f90](diffeq_singlestep.f90), which provides the adaptive step-size driver, dense output, and the common `solve` and `get_solution` interface.

The integrators are grouped below by source folder.

| Integrator | Type | Order | Embedded Estimate | Source |
|---|---|---|---|---|
| `runge_kutta_23` | Explicit | 3 | 2 | [explicit/diffeq_rk23.f90](explicit/diffeq_rk23.f90) |
| `runge_kutta_45` | Explicit | 5 | 4 | [explicit/diffeq_rk45.f90](explicit/diffeq_rk45.f90) |
| `tsitouras_54` | Explicit | 5 | 4 | [explicit/diffeq_tsitouras.f90](explicit/diffeq_tsitouras.f90) |
| `runge_kutta_853` | Explicit | 8 | 5 and 3 | [explicit/diffeq_rk853.f90](explicit/diffeq_rk853.f90) |
| `rosenbrock` | Linearly implicit | 4 | Yes | [implicit/diffeq_rosenbrock.f90](implicit/diffeq_rosenbrock.f90) |
| `kennedy_carpenter_4` | Implicit (ESDIRK) | 4 | 3 | [implicit/diffeq_kennedy_carpenter.f90](implicit/diffeq_kennedy_carpenter.f90) |
| `kennedy_carpenter_5` | Implicit (ESDIRK) | 5 | 4 | [implicit/diffeq_kennedy_carpenter.f90](implicit/diffeq_kennedy_carpenter.f90) |

### Explicit Methods

#### `runge_kutta_23`
The Bogacki--Shampine RK23 method. A third-order solution is paired with an embedded second-order solution, giving the local error estimate $e_{n+1} = y_{n+1}^{[3]} - y_{n+1}^{[2]}$. The method has four stages with first-same-as-last (FSAL) reuse. Its lower order reduces the cost per step, which makes it useful when high accuracy is not required.

#### `runge_kutta_45`
The Dormand--Prince RK45 method. A fifth-order solution is paired with an embedded fourth-order solution, and the fifth-order result is used after an accepted step. The local error estimate is $e_{n+1} = y_{n+1}^{[5]} - y_{n+1}^{[4]}$. The method has seven stages and FSAL reuse, so the final stage derivative becomes the first derivative of the next step. It is a general-purpose choice for moderate-to-high accuracy on non-stiff problems.

#### `tsitouras_54`
The Tsitouras fifth-order Runge--Kutta method with an embedded fourth-order estimator. The seven-stage FSAL tableau is designed for efficient adaptive integration of non-stiff systems.

#### `runge_kutta_853`
The Fehlberg RK8(5,3) method. It computes an eighth-order solution together with fifth- and third-order embedded solutions. The error estimate is $e_{n+1} = y_{n+1}^{[5]} - y_{n+1}^{[3]}$, while the accepted solution is the eighth-order result. The higher order is useful for smooth, non-stiff problems over long integration intervals.

### Implicit Methods

#### `rosenbrock`
A fourth-order Rosenbrock method for stiff systems. Rather than solving a nonlinear equation at every stage, it reuses one factored linear system based on the Jacobian $J = \partial f / \partial y$. Without a mass matrix the factored matrix is $A = \frac{1}{\gamma h}I - J$, and a user-supplied mass matrix produces $A = \frac{1}{\gamma h}M - J$. The mass matrix may be constant or state-dependent, in which case it is recomputed as needed. The method can solve index-1 DAEs when the initial conditions satisfy the algebraic constraints. Step sizes are chosen by a Rosenbrock-specific controller rather than the inherited PI controller. The iteration matrix is solved with pivoted LU, with a rank-revealing QR fallback for singular or unreliable systems.

#### `kennedy_carpenter_4` and `kennedy_carpenter_5`
Kennedy--Carpenter ESDIRK methods for stiff systems. `kennedy_carpenter_4` is the ARK4(3)6L[2]SA method, which is fourth order with a third-order embedded estimate. `kennedy_carpenter_5` is the ARK5(4)8L[2]SA method, which is fifth order with a fourth-order embedded estimate. Both share an implementation that solves each stage with a simplified Newton iteration. A supplied mass matrix is incorporated in each stage system as $M - h\gamma J$.

Every implicit stage shares the diagonal coefficient $\gamma$, so a single factorization serves all stages and Newton iterations within a step. The factorization is retained until the step size changes, and the Jacobian is retained across steps until the Newton iteration converges slowly or fails. A failed iteration first refreshes any stale Jacobian or factorization. If the iteration still fails, the step is rejected and retried with a smaller step size.

The Newton iteration limit and convergence tolerance are adjustable through `set_newton_step_limit` and `set_newton_tolerance`. The applied tolerance never exceeds one percent of the smallest positive integration tolerance. These solvers use the inherited, configurable PI step-size controller. Singular mass matrices require consistent initial conditions and an index-1 system.

### Choosing an Integrator

- **Non-stiff problems:** start with `runge_kutta_45` or `tsitouras_54`. Use `runge_kutta_853` for smooth problems that need high accuracy over long intervals, and `runge_kutta_23` when modest accuracy is acceptable at low cost.
- **Stiff problems or mass-matrix DAEs:** use `rosenbrock` or one of the Kennedy--Carpenter methods. The explicit methods are not intended for stiff systems or mass-matrix DAEs.

For a fuller comparison, see the [Choosing an Integrator](../../README.md#choosing-an-integrator) section of the main README.

### Examples

The examples below are located in the top-level [examples/](../../examples) folder.

#### Explicit vs. Implicit: Van der Pol Equation
[rk853_example.f90](../../examples/rk853_example.f90) and [rosenbrock_example.f90](../../examples/rosenbrock_example.f90) solve the same Van der Pol problem ($\mu = 5$, $0 \le x \le 50$, $y_1(0) = 2$, $y_2(0) = 0$). The first uses the explicit `runge_kutta_853` method, and the second uses the implicit `rosenbrock` method. The results can be compared directly, since both should produce the same solution curve. The explicit method is well suited to this non-stiff case, while the Rosenbrock method is designed for problems where explicit methods would need prohibitively small steps.

Explicit method, `runge_kutta_853`:
```fortran
type(runge_kutta_853) :: integrator
type(ode_container) :: mdl
real(real64), allocatable :: sol(:,:)

mdl%fcn => vanderpol
call integrator%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0])
sol = integrator%get_solution()
```

![](../../images/rk853_example.png?raw=true)

Implicit method, `rosenbrock`:
```fortran
type(rosenbrock) :: integrator
type(ode_container) :: mdl
real(real64) :: mu
real(real64), allocatable :: sol(:,:)

mu = 5.0d0
mdl%fcn => vanderpol_args
call integrator%solve(mdl, [0.0d0, 5.0d1], [2.0d0, 0.0d0], args = mu)
sol = integrator%get_solution()
```

![](../../images/rosenbrock_example.png?raw=true)

#### Differential-Algebraic Equations: Constrained Pendulum
[dae_example.f90](../../examples/dae_example.f90) models a pendulum as a mass-matrix DAE and solves it with `rosenbrock`. The mass matrix is supplied through `mdl%mass_matrix`, and the initial conditions satisfy the algebraic constraint. The model parameters are passed through `args`.

```fortran
type(rosenbrock) :: integrator
type(ode_container) :: mdl
type(cartesian_pendulum_properties) :: args

mdl%fcn => cartesian_pendulum
mdl%mass_matrix => cartesian_pendulum_mass_matrix
call mdl%set_is_mass_matrix_dependent(.false.)

call integrator%solve( &
    mdl, &
    [0.0d0, 1.0d1], &
    [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0], &
    args = args &
)
sol = integrator%get_solution()
```

![](../../images/dae_results.png?raw=true)