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
handle, and an output `c_ode_event_statistics` object. For `N` coordinates and `M` constraints, pass
`nstates = 2 * (N + M)`.  The solver requires the position and velocity values specified in the initial state vector satisfy the constraints; however, the solver will determine proper initial values for the multiplier components.

The callback declarations describe the exact array lengths and storage. In
particular, the DAE force callback receives the coordinate and velocity state
but writes only the physical force vector; it must not write multiplier
components.

Zero-initialize `c_ggl_container` before assigning its callbacks and parameters.
To enable event detection, set `events` to a `c_ode_event_function` callback
and `event_equation_count` to a positive count. The callback receives the full
augmented state, including both multiplier blocks. Integration stops at the
earliest detected event in integration order and includes its time and state
in the results. The statistics report the zero-based event index, or `-1` if
no event is found. A `NULL` callback or a count less than one disables event
detection.

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
