## C API

The C API exposes the library's integrators through the functions, callback types, and structures declared in [diffeq.h](diffeq.h). The Fortran implementation is in the `diffeq_c_*.f90` files in this folder and uses `iso_c_binding`. A C program includes `diffeq.h`, links against the `diffeq` library, initializes a settings structure, supplies function-pointer callbacks, and receives an opaque results handle.

| Function | Method | Mass matrix | Jacobian | Events | Source |
|---|---|---|---|---|---|
| `c_runge_kutta_23` | Bogacki--Shampine 3(2) | No | No | `_events` variant | [diffeq_c_integrators.f90](diffeq_c_integrators.f90) |
| `c_runge_kutta_45` | Dormand--Prince 5(4) | No | No | `_events` variant | |
| `c_runge_kutta_853` | Hairer 8(5,3) | No | No | `_events` variant | |
| `c_tsitouras` | Tsitouras 5(4) | No | No | `_events` variant | |
| `c_adams` | Adams PECE, orders 1--12 | Nonsingular | No | `_events` variant | |
| `c_rosenbrock` | Rosenbrock, 4th order | Singular allowed | Optional | `_events` variant | |
| `c_kennedy_carpenter_4` / `_5` | ESDIRK, 4th/5th order | Singular allowed | Optional | `_events` variant | |
| `c_bdf` | BDF, orders 1--5 | Singular allowed | Optional | `_events` variant | |
| `c_dae_rosenbrock`, `c_dae_kennedy_carpenter_4` / `_5`, `c_dae_bdf` | Stabilized index-1 GGL | Physical mass | Optional | Via `c_ggl_container` | [diffeq_c_index_reduction.f90](diffeq_c_index_reduction.f90) |

The sparse Jacobian and sparse mass-matrix paths are currently available only through the Fortran API.

### Callbacks and Storage

All callbacks receive dimensions explicitly and must write every output element:

| Type | Signature | Purpose |
|---|---|---|
| `c_ode` | `void (int n, double x, const double *y, double *dydx)` | Right-hand side $y'=f(x,y)$ |
| `c_ode_jacobian` | `void (int m, int n, double x, const double *y, double *jac, int ldj)` | Jacobian $\partial f/\partial y$ |
| `c_ode_mass_matrix` | `void (int n, double x, const double *y, double *m, int ldm)` | Mass matrix $M(x,y)$ |
| `c_ode_event_function` | `void (int neqn, int nstates, double x, const double *y, double *f)` | Event functions whose zeros are located |
| `c_dae_algebraic_constraints` | `void (int nc, int nstates, double x, const double *q, double *g)` | Position constraints $g(q)$ |
| `c_dae_constraint_curvature` | `void (int nc, int nstates, double x, const double *q, const double *direction, const double *v, double *crv)` | Curvature term $DG(q)[\dot q]v$ |

Matrices are stored in **column-major** order: element $(i,j)$ of a matrix with leading dimension `ld` is at index `j * ld + i`, using zero-based indices. This applies to Jacobian, mass-matrix, and constraint-gradient callbacks and to the solution matrix returned by `c_get_ode_results`. Pass `NULL` for any optional Jacobian or mass-matrix callback. The solver then uses a finite-difference Jacobian or an identity mass matrix.

### Settings, Results, and Errors

`c_default_ode_settings` fills a `c_ode_settings` structure with default values. Adjust individual fields before the solve:

- `absolute_tolerance` and `relative_tolerance` apply to every equation.
- `minimum_step_size`, `maximum_step_size`, and `step_limit` bound the step-size control.
- `beta` is the PI step-size controller exponent; zero disables PI control.
- `mass_matrix_constant` allows the mass matrix to be cached when it does not depend on the state.
- `allow_overshoot` permits a step past the final point, followed by interpolation back to it.
- `newton_step_limit` and `newton_tolerance` govern the Newton iteration in BDF and Kennedy--Carpenter.

Event solvers also take a `c_ode_event_solver_settings` structure, initialized with `c_default_ode_event_solver_settings`.

Each solver writes a `c_ode_results` handle that owns the solution and must be released with `c_free_ode_results`. With two entries in the independent-variable array, the solution contains every accepted step; with more entries, dense output returns the solution at exactly those points. Use the handle with:

- `c_get_ode_results_size` to obtain the row and column counts of the solution matrix.
- `c_get_ode_results` to copy the solution into a caller-allocated column-major matrix. The first column contains the independent variable, followed by one column for each state.
- `c_write_ode_results` to write the solution to a delimited text file.
- `c_get_ode_results_statistics` to obtain step, function-evaluation, Jacobian, and factorization counts.

Functions that return `int` return `C_NO_ERROR` (zero) on success. Other codes are `C_NULL_POINTER_ERROR`, `C_MATRIX_ROW_SIZE_ERROR`, `C_MATRIX_COLUMN_SIZE_ERROR`, `C_MATRIX_STRIDE_ERROR`, and `C_INVALID_INPUT_ERROR`. Check the code before using any output.

### Examples

The examples below are located in the top-level [examples/](../../examples) folder and are built with `BUILD_DIFFEQ_EXAMPLES=ON`. Each links against `diffeq` and adds this folder to its include path. The plotting examples also compile [c_plot_helper.f90](../../examples/c_plot_helper.f90), a small fplot wrapper declared in [c_plot_helper.h](../../examples/c_plot_helper.h). The excerpts omit some error handling for brevity; the full sources check every returned code.

#### Van der Pol Equation

[c_vanderpol_example.c](../../examples/c_vanderpol_example.c) solves the Van der Pol equation with $\mu = 5$ using `c_tsitouras`. Because only the start and end points are requested, every accepted step is returned. The example copies the solution into a matrix, writes it to a file, and plots $dy/dx$.

```c
void vanderpol(int n, double x, const double *y, double *dydx)
{
    const double mu = 5.0;
    dydx[0] = y[1];
    dydx[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
}

double x[] = {0.0, 50.0};
double ic[] = {2.0, 0.0};
int m, n, flag;
double *sol;
c_ode_results rst;
c_ode_settings settings;

c_default_ode_settings(&settings);
flag = c_tsitouras(2, 2, vanderpol, x, ic, &settings, &rst);
if (flag != C_NO_ERROR) return flag;

c_get_ode_results_size(rst, &m, &n);
sol = (double*)malloc((size_t)(m * n * sizeof(double)));
flag = c_get_ode_results(rst, m, n, sol, m);
flag = c_write_ode_results("example file.txt", rst, ',');

c_plot_xy(m, sol, &sol[2 * m], "x", "dy/dx");

free(sol);
c_free_ode_results(rst);
```

`&sol[2 * m]` points to the third column, which holds $dy/dx$.

![Van der Pol dy/dx from the C API](../../images/c_vanderpol_example.png?raw=true)

#### Dense Output

[c_dense_output_example.c](../../examples/c_dense_output_example.c) solves the same problem but requests 100 evenly spaced output points. The integrator still chooses its own adaptive steps and uses dense output to interpolate the solution at the requested points.

```c
const double xmax = 50.0;
const int npts = 100;
double x[npts];
double ic[] = {2.0, 0.0};

dx = xmax / (npts - 1.0);
for (i = 0; i < npts; ++i) x[i] = i * dx;

c_default_ode_settings(&settings);
flag = c_tsitouras(npts, 2, vanderpol, x, ic, &settings, &rst);
```

The result has 100 rows, independent of the number of internal steps. The example plots $y$ at these points.

![Van der Pol solution on a requested output grid](../../images/c_dense_output_example.png?raw=true)

#### Comparing Integrators

[c_comparison_example.c](../../examples/c_comparison_example.c) solves the Van der Pol problem with each integrator and prints its statistics. The implicit solvers take `NULL` for the Jacobian and mass-matrix callbacks, and setting `settings.beta` enables PI step-size control for an additional Runge--Kutta 4/5 solve.

```c
c_ode_results rst;
c_ode_results_statistics stats;

c_default_ode_settings(&settings);

flag = c_runge_kutta_45(2, 2, vanderpol, x, ic, &settings, &rst);
flag = c_get_ode_results_statistics(rst, &stats);
c_free_ode_results(rst);

flag = c_rosenbrock(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
flag = c_get_ode_results_statistics(rst, &stats);
c_free_ode_results(rst);

flag = c_bdf(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
flag = c_get_ode_results_statistics(rst, &stats);
c_free_ode_results(rst);

flag = c_adams(2, 2, vanderpol, x, ic, &settings, NULL, &rst);
flag = c_get_ode_results_statistics(rst, &stats);
c_free_ode_results(rst);

settings.beta = 0.1;
flag = c_runge_kutta_45(2, 2, vanderpol, x, ic, &settings, &rst);
settings.beta = 0.0;
```

This example has no plot. With default settings, it reports:

| Integrator | Successful steps | Failed steps | Function evaluations | Jacobian evaluations | Factorizations |
|---|---:|---:|---:|---:|---:|
| Runge--Kutta 2/3 | 2464 | 26 | 7471 | -- | -- |
| Runge--Kutta 4/5 | 582 | 87 | 4015 | -- | -- |
| Runge--Kutta 4/5 with PI control | 1106 | 0 | 6637 | -- | -- |
| Runge--Kutta 8/5/3 | 924 | 84 | 12097 | -- | -- |
| Tsitouras 4/5 | 521 | 102 | 4362 | -- | -- |
| Rosenbrock | 1100 | 37 | 11260 | 1100 | 1137 |
| Kennedy--Carpenter 4 | 482 | 145 | 11263 | 207 | 559 |
| Kennedy--Carpenter 5 | 245 | 73 | 9278 | 173 | 315 |
| BDF | 1583 | 59 | 16331 | 3276 | 3276 |
| Adams | 1132 | 70 | 2336 | -- | -- |

These counters measure solver work, not elapsed time.

#### Event Detection

[c_event_detection_example.c](../../examples/c_event_detection_example.c) models a ball dropped from 5 m. The event function returns the height, so the integration stops at ground contact. A second solve restarts from the impact time with zero height and 70% of the impact speed, directed upward.

```c
void ball(int n, double t, const double *y, double *dydt)
{
    const double g = 9.81;
    dydt[0] = y[1];
    dydt[1] = -g;
}

void ball_ground_constraint(int neqn, int nstates, double t, const double *y, double *f)
{
    f[0] = y[0];
}

double t[] = {0.0, 10.0};
double ic[] = {5.0, 0.0};
c_ode_event_solver_settings esettings;
c_ode_event_statistics e1, e2;

c_default_ode_settings(&settings);
c_default_ode_event_solver_settings(&esettings);
settings.maximum_step_size = 1.0e-2;

flag = c_tsitouras_events(2, 2, 1, ball, t, ic, &settings,
    &esettings, ball_ground_constraint, &r1, &e1);

/* Copy r1 into s1 with c_get_ode_results, then restart from the impact. */
t[0] = s1[INDEX(m1 - 1, 0, m1)];
ic[0] = 0.0;
ic[1] = -k * s1[INDEX(m1 - 1, 2, m1)];
flag = c_tsitouras_events(2, 2, 1, ball, t, ic, &settings,
    &esettings, ball_ground_constraint, &r2, &e2);
```

`INDEX(i, j, ld)` expands to `j * ld + i`, and `k = 0.7` is the coefficient of restitution. The `c_ode_event_statistics` structure reports the zero-based index of the triggered event, or `-1` if none occurred, along with event-function evaluations and root-finding iterations. A `NULL` event function or an event count below one disables event detection.

The maximum step size only produces a smoother plot:

![Bouncing ball with a 0.01 s maximum step](../../images/c_event_detection_example.png?raw=true)

With default settings, the integrator takes large steps on this simple problem and still locates the events accurately:

![Bouncing ball with default settings](../../images/c_event_detection_example_default.png?raw=true)

#### Stabilized Pendulum

[c_stabilized_pendulum_example.c](../../examples/c_stabilized_pendulum_example.c) solves a Cartesian pendulum with `c_ggl_container` and `c_dae_bdf`. The container describes a holonomically constrained mechanical model, which the library reduces to a stabilized index-1 DAE. See the [DAE index reduction README](../dae/README.md) for the formulation.

```c
const double body_mass = 2.0;
const double length = 1.5;

void force(int n, double t, const double *x, double *frc)
{
    frc[0] = 0.0;
    frc[1] = -body_mass * 9.81;
}

void physical_mass(int n, double t, const double *x, double *m, int ldm)
{
    m[0] = body_mass;
    m[1] = 0.0;
    m[2] = 0.0;
    m[3] = body_mass;
}

void constraints(int nc, int nstates, double t, const double *x, double *g)
{
    g[0] = x[0] * x[0] + x[1] * x[1] - length * length;
}

void gradient(int m, int n, double t, const double *x, double *grad, int ldg)
{
    int i;
    for (i = 0; i < n; ++i) grad[i] = 2.0 * x[i];
}

void curvature(int nc, int nstates, double t, const double *q,
    const double *direction, const double *v, double *crv)
{
    crv[0] = 2.0 * (direction[0] * v[0] + direction[1] * v[1]);
}

c_ggl_container model = {0};
double ic[] = {length, 0.0, 0.0, 0.0, 0.0, 0.0};
double t[] = {0.0, 10.0};
c_ode_event_statistics estats;

model.alpha = 2.0;
model.beta = 3.0;
model.constraint_count = 1;
model.constraints = constraints;
model.curvature = curvature;
model.force = force;
model.gradient = gradient;
model.mass = physical_mass;
model.jacobian = NULL;

c_default_ode_settings(&settings);
flag = c_dae_bdf(&model, 2, 6, t, ic, &settings, &rst, &estats);
```

For $N$ coordinates and $M$ constraints, `nstates` is $2(N+M)$ and the state is ordered as $[q, v, \lambda, \mu]$; here it is `[qx, qy, vx, vy, lambda, mu]`. The `force` callback receives `[q, v]` but writes only the $N$ applied forces. Zero-initializing the container leaves the optional `jacobian` and `events` callbacks disabled.

Supply physically consistent positions and velocities. The wrapper computes consistent multipliers and overwrites those entries of `ic`, but it does not correct $q$ or $v$. `model.alpha` and `model.beta` damp position and velocity constraint drift; they are unrelated to the PI-controller setting `settings.beta`. The same container can be passed to `c_dae_rosenbrock`, `c_dae_kennedy_carpenter_4`, or `c_dae_kennedy_carpenter_5`.

![Stabilized pendulum coordinates from the C API](../../images/c_stabilized_pendulum_example.png?raw=true)
