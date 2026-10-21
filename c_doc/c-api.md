# DIFFEQ C API

The DIFFEQ C API provides C-compatible access to the library's ordinary
differential equation (ODE) solvers and its stabilized index-1 DAE solvers.
The public declarations are in `diffeq.h`; the generated reference documents
each callback, data structure, and function.

## Getting started

Include `diffeq.h` and link your program to the installed DIFFEQ library. For a
CMake consumer, link the exported target:

```cmake
find_package(diffeq CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE diffeq::diffeq)
```

Initialize settings with `c_default_ode_settings` before changing any fields.
Solver functions return a `c_ode_results` handle on success. Query its size,
copy its solution, and release it with `c_free_ode_results`:

```c
c_ode_settings settings;
c_ode_results results;
int rows, columns;

c_default_ode_settings(&settings);
if (c_runge_kutta_45(point_count, state_count, rhs, times, initial_state,
        &settings, &results) != C_NO_ERROR) {
    /* Handle the integration error. */
}

if (c_get_ode_results_size(results, &rows, &columns) == C_NO_ERROR) {
    /* Allocate rows * columns doubles, then copy with leading dimension rows. */
}
c_free_ode_results(results);
```

Use the function declaration reference for exact argument requirements and
error codes. Output handles must be freed even when the caller no longer needs
the computed solution.

## Conventions

- Solution matrices and callback Jacobians use column-major storage. For a
  matrix with leading dimension `ld`, element `(row, column)` is stored at
  `column * ld + row`, using zero-based C indices.
- Callback array lengths and matrix dimensions are provided explicitly.
  Callbacks must write every output element required by their declaration.
- Solver time grids contain at least two points. Requested points are returned
  in order; with more than two points the solver will use dense output.
- Optional ODE Jacobian and mass-matrix callbacks may be `NULL` where allowed
  by the solver declaration. The solvers use finite-difference Jacobians when no analytic
  Jacobian is supplied.
- Utilize a `c_ode_results` handle to retrieve the integration results, and release the handle after finished with the results by calling `c_free_ode_results`.

## Differential-algebraic equations

The C DAE solvers use the stabilized Gear-Gupta-Leimkuhler (GGL) formulation
for holonomically constrained mechanical systems. A @ref c_ggl_container
describes the applied force, physical
mass matrix, constraints, constraint gradient, curvature, and stabilization
parameters. For `N` coordinates and `M` constraints, the integrated state has
`2 * (N + M)` entries ordered as coordinates, velocities, position
multipliers, and velocity multipliers.

The `constraints` callback returns the position residuals. The `gradient`
callback returns the constraint Jacobian in column-major order, and
`curvature` returns the directional curvature terms without the acceleration
term. These callbacks must accept perturbed, off-constraint states. The
physical mass callback may be omitted to use the identity matrix; set
`settings.mass_matrix_constant` to false when a supplied mass matrix varies
with time or coordinates.

Use @ref c_dae_rosenbrock, @ref c_dae_kennedy_carpenter_4,
@ref c_dae_kennedy_carpenter_5, or @ref c_dae_bdf to integrate a GGL model.
These methods require a consistent initial state, including the multiplier
entries. The complete model and callback contracts are documented with
@ref c_ggl_container in the generated reference.

For shared details about solver arguments, result layout, and error handling,
see @ref c_api_reference "the C API reference notes".

## Reference

The generated reference includes all ODE and DAE solver functions, result
accessors, settings, callback types, and model structures declared in
`diffeq.h`.
