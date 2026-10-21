@page c_api_reference DIFFEQ C API Reference Notes

The generated declaration reference for `diffeq.h` is the authoritative source
for C types, parameter lists, callback signatures, and return codes. This page
summarizes conventions shared by the routines.

## Integrator arguments

ODE integrators take the number of requested points, state count, an ODE
callback, the requested time grid, initial state, settings, and an output
results handle. Rosenbrock, Kennedy-Carpenter, and BDF ODE interfaces also
accept optional Jacobian and mass-matrix callbacks.

DAE integrators instead take a `c_ggl_container`, requested point count, full
augmented-state count, time grid, mutable initial state, settings, and output
handle. For `N` coordinates and `M` constraints, pass
`nstates = 2 * (N + M)`.  The solver requires the position and velocity values specified in the initial state vector satisfy the constraints; however, the solver will determine proper initial values for the multiplier components.

The callback declarations describe the exact array lengths and storage. In
particular, the DAE force callback receives the coordinate and velocity state
but writes only the physical force vector; it must not write multiplier
components.

## Results and errors

`c_get_ode_results_size` reports the dimensions required by
`c_get_ode_results`. The output matrix has one column for time followed by one
column per state component, and is stored column-major. Its leading dimension
must be at least the number of rows. `c_get_ode_results_statistics` returns
step, function, Jacobian, and factorization counts.

Functions returning `C_NO_ERROR` succeeded. Other return values are defined by
the error constants in `diffeq.h`. Check the status before reading an output
handle or result buffer, and release a successful result handle with
`c_free_ode_results`.
