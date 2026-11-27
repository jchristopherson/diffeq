## DAE Index Reduction

Mechanical systems with holonomic constraints $g(q)=0$ lead to index-3 differential-algebraic equations (DAEs), which the integrators in this library cannot solve directly. The `ggl_container` type in [diffeq_index_reduction.f90](diffeq_index_reduction.f90) reduces such a model to a stabilized index-1 DAE based on the Gear--Gupta--Leimkuhler (GGL) formulation. It extends `ode_container`, so the reduced system can be passed to any integrator that accepts a singular mass matrix.

| Integrator | Dense path | Sparse path | Source |
|---|---|---|---|
| `bdf` | Yes | Yes | [diffeq_bdf.f90](../multistep/diffeq_bdf.f90) |
| `rosenbrock` | Yes | No | [diffeq_rosenbrock.f90](../singlestep/implicit/diffeq_rosenbrock.f90) |
| Kennedy--Carpenter ESDIRK (4th and 5th order) | Yes | No | [diffeq_kennedy_carpenter.f90](../singlestep/implicit/diffeq_kennedy_carpenter.f90) |

### Stabilized Index-1 GGL Formulation

For $N$ coordinates $q$, physical velocities $v$, and $M$ independent constraints, the container adds two multiplier vectors: $\lambda$, which enforces the velocity-level constraint through the constraint forces, and $\mu$, which enforces the position constraint by correcting the kinematic relation. With $G(q)=\partial g/\partial q$, the augmented state and mass matrix are

$$
z=[q_N,v_N,\lambda_M,\mu_M]^T, \qquad
\mathcal M=\operatorname{diag}(I_N,M_N,0_M,0_M),
$$

and the container supplies the equations

$$
\begin{aligned}
\dot q &= v+G^T\mu,\\
M\dot v &= F(x,q,v)-G^T\lambda,\\
0 &= G(v+G^T\mu)+\alpha g(q),\\
0 &= DG(q)[v+G^T\mu]\,v+GM^{-1}(F-G^T\lambda)+\beta Gv.
\end{aligned}
$$

The last two rows replace the original position and velocity constraints with their stabilized evolution equations, so that

$$
\dot g=-\alpha g, \qquad \frac{d(Gv)}{dx}=-\beta Gv.
$$

Both constraints remain zero in the exact solution when they are initially satisfied. In that case $\mu=0$, and $\lambda$ equals the classical constraint-force multiplier. Positive gains $\alpha$ and $\beta$ damp numerical drift; zero gains preserve, rather than correct, any initial residual. Gains default to one and have inverse-time units when $x$ is time. Very large gains can add stiffness. Stabilization is not projection: accepted steps and interpolated output are not forced onto the constraint surface.

This is a stabilized index-1 variant inspired by GGL, not the classical index-2 GGL system. The augmented mass matrix is singular, but the physical mass $M$ must be nonsingular. Independent constraints and a symmetric positive-definite $M$ make $GG^T$ and $GM^{-1}G^T$ nonsingular. Explicitly time-dependent constraints are not supported.

### Defining a Model

The callbacks have shapes specific to the mechanical model, not to the augmented state:

| Callback | Input | Output |
|---|---|---|
| `fcn` (required) | $[q,v]$, length $2N$ | Applied force $F$, length $N$, excluding constraint forces |
| `constraints` (required) | $q$, length $N$ | Position residual $g$, length $M$ |
| `constraint_jacobian` (required) | $q$, length $N$ | $G$, an $M\times N$ matrix |
| `constraint_curvature` (required) | $q$, direction $\dot q$, and velocity $v$ | $DG(q)[\dot q]v$, length $M$; zeros for linear constraints |
| `mass_matrix` (optional) | $q$, length $N$ | Dense physical $N\times N$ mass |
| `sparse_mass_matrix` (optional) | $q$, length $N$ | Physical $N\times N$ mass in CSR format |
| `jacobian` / `sparse_jacobian` (optional) | Augmented state $z$ | Jacobian of the complete augmented RHS |

The constraint callbacks must accept off-manifold coordinates, because the solvers evaluate perturbed and predicted states. There is no finite-difference fallback for $G$ or the curvature term. Call `set_is_mass_matrix_dependent(.true.)` when the physical mass depends on time or coordinates. The inherited `fcn_ad` is not supported, because differentiating only the applied force would omit the constraint terms.

Configure the container with `initialize(N, M, alpha, beta)` or the `ggl_container(N, M, alpha, beta)` constructor; `set_stabilization_parameters` changes the gains later. Before solving:

1. `validate_initial_conditions(x, z, tolerance)` checks that $g(q)$ and $Gv$ are within an absolute tolerance. It does not change the state.
2. `initialize_multipliers(x, z)` computes consistent $\lambda$ and $\mu$ for the supplied $q$ and $v$. It does not correct inconsistent positions or velocities.

### Dense and Sparse Systems

The container always resolves the physical mass internally when it forms the acceleration $M^{-1}(F-G^T\lambda)$. The integrator, meanwhile, works with the augmented mass matrix and the Jacobian of the augmented RHS.

- **Dense systems:** the physical mass comes from `mass_matrix` or defaults to the identity. Physical-mass solves use LU factorization. If no augmented Jacobian is supplied, the inherited finite-difference implementation differentiates the complete augmented RHS.
- **Sparse systems:** a `sparse_mass_matrix` callback supplies the physical mass in CSR format. The container solves with it using ILUT-preconditioned GMRES and assembles the sparse augmented matrix $\operatorname{diag}(I,M,0,0)$ itself. If both mass callbacks are assigned, the dense callback takes precedence. Dense-only integrators still accept a sparse mass callback; it is converted to a dense matrix when they request one.

`bdf` selects its sparse Newton iteration when `has_sparse_jacobian()` returns true, either because `sparse_jacobian` is assigned or because a derived type overrides the method. The sparse Jacobian must differentiate the complete augmented RHS, not just the applied force; automatic sparse Jacobian construction is not yet available. Extending `ggl_container` and overriding `has_sparse_jacobian` and `compute_sparse_jacobian` gives the Jacobian access to the container's `compute_rhs`, as the sparse example below demonstrates.

Factorization and RHS scratch storage are retained inside the container and reused across calls. A container must not be used by concurrent or nested integrations.

### Examples

The examples below are located in the top-level [examples/](../../examples) folder and use fplot to display their results.

#### Dense GGL: Stabilized Pendulum

[stabilized_pendulum_example.f90](../../examples/stabilized_pendulum_example.f90) models a planar pendulum in Cartesian coordinates with a single length constraint, $g(q)=q^Tq-L^2$. Its gradient is $G=2q^T$, and its curvature term is $DG(q)[\dot q]v=2\dot q^Tv$. The model uses a dense physical mass matrix and the dense `bdf` path. The state ordering is $[q_x,q_y,v_x,v_y,\lambda,\mu]$.

```fortran
type(ggl_container) :: model
type(bdf) :: solver
real(real64), dimension(6) :: initial_state
real(real64), allocatable, dimension(:,:) :: sol

model = ggl_container(2, 1, alpha = 2.0d0, beta = 3.0d0)
model%fcn => force
model%mass_matrix => physical_mass
model%constraints => position_constraint
model%constraint_jacobian => constraint_gradient
model%constraint_curvature => constraint_curvature
call model%set_is_mass_matrix_dependent(.false.)

initial_state = [length, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
if (.not.model%validate_initial_conditions(0.0d0, initial_state, 1.0d-10)) &
    error stop "Inconsistent physical initial conditions"
call model%initialize_multipliers(0.0d0, initial_state)

call solver%set_absolute_tolerance(1.0d-9)
call solver%set_relative_tolerance(1.0d-9)
call solver%solve(model, [0.0d0, 1.0d1], initial_state)
sol = solver%get_solution()
```

The callbacks are defined in the example. The figure shows the Cartesian coordinates and the corresponding pendulum angle.

![Stabilized pendulum coordinates and angle](../../images/stabilized_pendulum_example.png?raw=true)

#### Sparse GGL: Falling Chain

[sparse_ggl_example.f90](../../examples/sparse_ggl_example.f90) models a planar chain of 30 point masses joined by rigid, massless links, with the first link pinned at the origin. Each link contributes the constraint $g_k=\|p_k-p_{k-1}\|^2-L^2$. This gives 60 coordinates, 30 constraints, and a 180-equation augmented system. Because each constraint couples only neighbouring particles, fewer than 5% of the augmented Jacobian entries are nonzero.

The physical mass is supplied only as a diagonal CSR matrix. The `chain_model` type extends `ggl_container` and overrides the sparse-Jacobian methods, which selects the sparse `bdf` path:

```fortran
type, extends(ggl_container) :: chain_model
    integer :: jacobian_calls = 0
contains
    procedure :: has_sparse_jacobian => chain_has_sparse_jacobian
    procedure :: compute_sparse_jacobian => chain_sparse_jacobian
end type
```

The overridden `compute_sparse_jacobian` uses colored central differences of `compute_rhs`. Columns that cannot affect the same rows are perturbed together, so each Jacobian requires 60 RHS evaluations instead of 360. The model is then configured and solved like the dense case:

```fortran
type(chain_model) :: model
type(bdf) :: solver

call model%initialize(ncoord, nlinks, alpha = alpha, beta = beta)
model%fcn => chain_force
model%sparse_mass_matrix => chain_mass
model%constraints => chain_constraints
model%constraint_jacobian => chain_gradient
model%constraint_curvature => chain_curvature
call model%set_is_mass_matrix_dependent(.false.)

if (.not.model%validate_initial_conditions(0.0d0, state, 1.0d-12)) &
    error stop "Inconsistent initial conditions"
call model%initialize_multipliers(0.0d0, state)

call solver%set_absolute_tolerance(1.0d-7)
call solver%set_relative_tolerance(1.0d-7)
call solver%solve(model, times, state)
sol = solver%get_solution()
```

The chain starts straight and at rest, 60 degrees from vertical. Before solving, the example compares the sparse Jacobian with the inherited dense finite-difference Jacobian. Afterwards, it checks that the sparse Jacobian was used and reports the constraint drift and energy change. With tolerances of $10^{-7}$, the maximum position drift $\max|g(q)|$ is about $2\times10^{-9}$, and the maximum velocity drift $\max|G(q)v|$ is about $4\times10^{-7}$. The left plot shows snapshots of the falling chain; the right plot shows both drift measures on a logarithmic scale.

![Falling chain snapshots and constraint drift](../../images/sparse_ggl_example.png?raw=true)
