/*
* Copyright (c) 2023-2026 Jason Christopherson
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include "diffeq.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { coordinate_count = 2, constraint_count = 1, state_count = 6 };

typedef int (*dae_integrator)(int, int, const c_ggl_container *,
    const double *, double *, const c_ode_settings *, c_ode_results *,
    c_ode_event_statistics *);

static const char *test_name;
static int force_calls;
static int mass_calls;
static int constraint_calls;
static int gradient_calls;
static int curvature_calls;
static int event_calls;

static void check(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "C DAE integrator test failed [%s]: %s\n",
            test_name, message);
        exit(EXIT_FAILURE);
    }
}

static void check_close(double actual, double expected, double tolerance,
    const char *message)
{
    if (!(fabs(actual - expected) < tolerance)) {
        fprintf(stderr, "Expected %.17g, received %.17g\n", expected, actual);
        check(0, message);
    }
}

static void oscillator_force(int n, double x, const double *y, double *force)
{
    (void)x;
    check(n == 2 * coordinate_count, "force state count");
    ++force_calls;
    force[0] = -y[0];
    force[1] = 0.0;
}

static void oscillator_events(int nevents, int nstates, double x,
    const double *y, double *values)
{
    check(nevents == 2 && nstates == state_count, "event callback dimensions");
    check(isfinite(x), "event callback time");
    ++event_calls;
    values[0] = y[0] - cos(0.8);
    values[1] = y[0] - cos(0.6);
}

static void noncrossing_events(int nevents, int nstates, double x,
    const double *y, double *values)
{
    oscillator_events(nevents, nstates, x, y, values);
    values[0] = y[0] + 2.0;
    values[1] = y[0] + 3.0;
}

static void identity_mass(int n, double x, const double *q, double *mass,
    int ldm)
{
    (void)x;
    (void)q;
    check(n == coordinate_count && ldm >= n, "mass matrix dimensions");
    ++mass_calls;
    mass[0] = 1.0;
    mass[1] = 0.0;
    mass[ldm] = 0.0;
    mass[ldm + 1] = 1.0;
}

static void constrain_y(int nconstraints, int nstates, double x,
    const double *q, double *constraints)
{
    (void)x;
    check(nconstraints == constraint_count && nstates == coordinate_count,
        "constraint dimensions");
    ++constraint_calls;
    constraints[0] = q[1];
}

static void constrain_y_gradient(int m, int n, double x, const double *q,
    double *gradient, int ldj)
{
    (void)x;
    (void)q;
    check(m == constraint_count && n == coordinate_count && ldj >= m,
        "constraint gradient dimensions");
    ++gradient_calls;
    gradient[0] = 0.0;
    gradient[ldj] = 1.0;
}

static void linear_constraint_curvature(int nconstraints, int nstates,
    double x, const double *q, const double *direction, const double *velocity,
    double *curvature)
{
    (void)x;
    (void)q;
    (void)direction;
    (void)velocity;
    check(nconstraints == constraint_count && nstates == coordinate_count,
        "constraint curvature dimensions");
    ++curvature_calls;
    curvature[0] = 0.0;
}

static int dae_rosenbrock(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output, c_ode_event_statistics *estats)
{
    return c_dae_rosenbrock(model, nx, nstates, x, initial, settings, output,
        estats);
}

static int dae_kc4(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output, c_ode_event_statistics *estats)
{
    return c_dae_kennedy_carpenter_4(model, nx, nstates, x, initial, settings,
        output, estats);
}

static int dae_kc5(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output, c_ode_event_statistics *estats)
{
    return c_dae_kennedy_carpenter_5(model, nx, nstates, x, initial, settings,
        output, estats);
}

static int dae_bdf(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output, c_ode_event_statistics *estats)
{
    return c_dae_bdf(model, nx, nstates, x, initial, settings, output, estats);
}

static const struct {
    const char *name;
    dae_integrator solve;
} integrators[] = {
    {"dae_rosenbrock", dae_rosenbrock},
    {"dae_kc4", dae_kc4},
    {"dae_kc5", dae_kc5},
    {"dae_bdf", dae_bdf}
};

static void test_solution(dae_integrator solve, c_ode_event_function events,
    int nevents)
{
    const double times[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    double initial[state_count] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double matrix[5 * (state_count + 1)];
    c_ggl_container model = {0};
    c_ode_settings settings;
    c_ode_results results = NULL;
    c_ode_results_statistics stats;
    c_ode_event_statistics estats = {-1, -1, -2};
    int rows, cols, i;

    model.force = oscillator_force;
    model.mass = identity_mass;
    model.constraints = constrain_y;
    model.gradient = constrain_y_gradient;
    model.curvature = linear_constraint_curvature;
    model.alpha = 2.0;
    model.beta = 3.0;
    model.constraint_count = constraint_count;
    model.events = events;
    model.event_equation_count = nevents;

    c_default_ode_settings(&settings);
    settings.absolute_tolerance = 1.0e-9;
    settings.relative_tolerance = 1.0e-9;
    settings.maximum_step_size = 0.05;
    force_calls = mass_calls = constraint_calls = gradient_calls =
        curvature_calls = 0;
    event_calls = 0;

    check(solve(5, state_count, &model, times, initial, &settings, &results,
        &estats) ==
        C_NO_ERROR, "solve constrained oscillator");
    check(results != NULL, "integration returns a result handle");
    check(c_get_ode_results_size(results, &rows, &cols) == C_NO_ERROR,
        "get result dimensions");
    check(rows == 5 && cols == state_count + 1, "result dimensions");
    check(c_get_ode_results(results, rows, cols, matrix, rows) == C_NO_ERROR,
        "copy result matrix");

    for (i = 0; i < rows; ++i) {
        double t = times[i];
        check_close(matrix[i], t, 1.0e-12, "reported time");
        check_close(matrix[i + rows], cos(t), 2.0e-6,
            "constrained oscillator position");
        check_close(matrix[i + 2 * rows], 0.0, 2.0e-7,
            "algebraic position constraint");
        check_close(matrix[i + 3 * rows], -sin(t), 2.0e-6,
            "constrained oscillator velocity");
        check_close(matrix[i + 4 * rows], 0.0, 2.0e-7,
            "algebraic velocity constraint");
    }

    check(c_get_ode_results_statistics(results, &stats) == C_NO_ERROR,
        "get integration statistics");
    check(stats.successful_steps > 0 && stats.function_evaluations > 0,
        "integration statistics are populated");
    check(force_calls > 0 && mass_calls > 0, "force and mass callbacks invoked");
    check(constraint_calls > 0 && gradient_calls > 0 && curvature_calls > 0,
        "constraint callbacks invoked");
    check(events && nevents > 0 ? event_calls > 0 : event_calls == 0,
        "event callback invoked only when enabled");
    check(estats.event_index == -1 && estats.iteration_count == 0 &&
        estats.function_evaluations == 0, "no-event statistics");
    c_free_ode_results(results);
}

static void test_events(dae_integrator solve, int reverse, int dense)
{
    double times[5];
    double initial[state_count] = {0};
    double *matrix;
    c_ggl_container model = {0};
    c_ode_settings settings;
    c_ode_event_statistics estats = {-1, -1, -2};
    c_ode_results results = NULL;
    const double root = reverse ? 0.8 : 0.6;
    int rows, cols, i;
    int nx = dense ? 5 : 2;

    for (i = 0; i < nx; ++i) {
        times[i] = reverse ? 1.0 - (double)i / (nx - 1) :
            (double)i / (nx - 1);
    }
    initial[0] = cos(times[0]);
    initial[2] = -sin(times[0]);
    model.force = oscillator_force;
    model.mass = identity_mass;
    model.constraints = constrain_y;
    model.gradient = constrain_y_gradient;
    model.curvature = linear_constraint_curvature;
    model.alpha = 2.0;
    model.beta = 3.0;
    model.constraint_count = constraint_count;
    model.events = oscillator_events;
    model.event_equation_count = 2;
    c_default_ode_settings(&settings);
    settings.absolute_tolerance = 1.0e-9;
    settings.relative_tolerance = 1.0e-9;
    settings.maximum_step_size = 0.05;
    event_calls = 0;

    check(solve(nx, state_count, &model, times, initial, &settings, &results,
        &estats) == C_NO_ERROR, "solve with events");
    check(results != NULL, "event integration returns a result handle");
    check(event_calls > 0, "event callback invoked");
    check(estats.event_index == (reverse ? 0 : 1),
        "zero-based earliest event index in integration order");
    check(estats.iteration_count > 0 && estats.function_evaluations > 0,
        "event statistics populated");
    check(c_get_ode_results_size(results, &rows, &cols) == C_NO_ERROR,
        "get event result dimensions");
    check(rows >= 2 && cols == state_count + 1, "event result dimensions");
    if (dense) check(rows == (reverse ? 2 : 4),
        "requested points before event are preserved");
    matrix = malloc((size_t)rows * cols * sizeof(*matrix));
    check(matrix != NULL, "allocate event solution");
    check(c_get_ode_results(results, rows, cols, matrix, rows) == C_NO_ERROR,
        "copy event solution");
    check_close(matrix[rows - 1], root, 2.0e-6, "event time");
    for (i = 0; i < rows; ++i) {
        double t = matrix[i];
        check(reverse ? t >= root - 2.0e-6 : t <= root + 2.0e-6,
            "output stops at earliest event");
        if (i > 0) check(reverse ? t < matrix[i - 1] : t > matrix[i - 1],
            "event output ordered in integration direction");
        if (dense && i < rows - 1) check_close(t, times[i], 1.0e-12,
            "pre-event requested time");
        check_close(matrix[i + rows], cos(t), 2.0e-6, "event position");
        check_close(matrix[i + 3 * rows], -sin(t), 2.0e-6, "event velocity");
        check_close(matrix[i + 2 * rows], 0.0, 2.0e-7,
            "event position constraint");
        check_close(matrix[i + 4 * rows], 0.0, 2.0e-7,
            "event velocity constraint");
        check_close(matrix[i + 5 * rows], 0.0, 2.0e-7,
            "event position multiplier");
        check_close(matrix[i + 6 * rows], 0.0, 2.0e-7,
            "event velocity multiplier");
    }
    free(matrix);
    c_free_ode_results(results);
}

int main(int argc, char **argv)
{
    size_t i;
    int matched = 0;

    check(argc <= 2, "expected at most one DAE integrator name");
    for (i = 0; i < sizeof(integrators) / sizeof(integrators[0]); ++i) {
        if (argc == 1 || strcmp(argv[1], integrators[i].name) == 0) {
            matched = 1;
            test_name = integrators[i].name;
            test_solution(integrators[i].solve, NULL, 0);
            test_solution(integrators[i].solve, NULL, 2);
            test_solution(integrators[i].solve, oscillator_events, 0);
            test_solution(integrators[i].solve, oscillator_events, -1);
            test_solution(integrators[i].solve, noncrossing_events, 2);
            test_events(integrators[i].solve, 0, 0);
            test_events(integrators[i].solve, 0, 1);
            test_events(integrators[i].solve, 1, 0);
            test_events(integrators[i].solve, 1, 1);
        }
    }
    test_name = "arguments";
    check(matched, "unknown DAE integrator name");
    return EXIT_SUCCESS;
}
