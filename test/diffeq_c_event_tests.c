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

typedef int (*event_integrator_function)(int, int, int, c_ode, const double *,
    const double *, const c_ode_settings *,
    const c_ode_event_solver_settings *, c_ode_event_function,
    c_ode_results *, c_ode_event_statistics *);

static const char *test_context;
static int callback_errors;
static int event_callback_calls;

static void check(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "C event test failed [%s]: %s\n", test_context, message);
        exit(EXIT_FAILURE);
    }
}

static void check_close(double actual, double expected, double tolerance,
    const char *message)
{
    if (!(fabs(actual - expected) <= tolerance)) {
        fprintf(stderr, "Expected %.17g, received %.17g\n", expected, actual);
        check(0, message);
    }
}

static void unit_rate(int n, double x, const double *y, double *dydx)
{
    int i;
    (void)x;
    (void)y;
    for (i = 0; i < n; ++i) dydx[i] = 1.0;
}

static void two_events(int neqn, int nstates, double x, const double *y,
    double *values)
{
    (void)x;
    ++event_callback_calls;
    if (neqn != 2 || nstates != 1) {
        callback_errors = 1;
        return;
    }
    values[0] = y[0] - 0.7;
    values[1] = y[0] - 0.3;
}

#define DEFINE_EVENT_ADAPTER(name, function) \
    static int name(int nx, int nstates, int nevents, c_ode fcn, \
        const double *x, const double *initial, const c_ode_settings *settings, \
        const c_ode_event_solver_settings *esettings, c_ode_event_function efcn, \
        c_ode_results *output, c_ode_event_statistics *estats) \
    { \
        return function(nx, nstates, nevents, fcn, x, initial, settings, \
            esettings, efcn, output, estats); \
    }

#define DEFINE_STIFF_EVENT_ADAPTER(name, function) \
    static int name(int nx, int nstates, int nevents, c_ode fcn, \
        const double *x, const double *initial, const c_ode_settings *settings, \
        const c_ode_event_solver_settings *esettings, c_ode_event_function efcn, \
        c_ode_results *output, c_ode_event_statistics *estats) \
    { \
        return function(nx, nstates, nevents, fcn, x, initial, settings, \
            esettings, efcn, NULL, NULL, output, estats); \
    }

DEFINE_EVENT_ADAPTER(rk23_events, c_runge_kutta_23_events)
DEFINE_EVENT_ADAPTER(rk45_events, c_runge_kutta_45_events)
DEFINE_EVENT_ADAPTER(rk853_events, c_runge_kutta_853_events)
DEFINE_EVENT_ADAPTER(tsitouras_events, c_tsitouras_events)
DEFINE_STIFF_EVENT_ADAPTER(rosenbrock_events, c_rosenbrock_events)
DEFINE_STIFF_EVENT_ADAPTER(kc4_events, c_kennedy_carpenter_4_events)
DEFINE_STIFF_EVENT_ADAPTER(kc5_events, c_kennedy_carpenter_5_events)
DEFINE_STIFF_EVENT_ADAPTER(bdf_events, c_bdf_events)

static int adams_events(int nx, int nstates, int nevents, c_ode fcn,
    const double *x, const double *initial, const c_ode_settings *settings,
    const c_ode_event_solver_settings *esettings, c_ode_event_function efcn,
    c_ode_results *output, c_ode_event_statistics *estats)
{
    return c_adams_events(nx, nstates, nevents, fcn, x, initial, settings,
        esettings, efcn, NULL, output, estats);
}

static const struct {
    const char *name;
    event_integrator_function solve;
} integrators[] = {
    {"rk23", rk23_events},
    {"rk45", rk45_events},
    {"rk853", rk853_events},
    {"tsitouras", tsitouras_events},
    {"adams", adams_events},
    {"rosenbrock", rosenbrock_events},
    {"kc4", kc4_events},
    {"kc5", kc5_events},
    {"bdf", bdf_events}
};

static void test_event_detection(event_integrator_function solve)
{
    const double x[] = {0.0, 1.0};
    const double initial[] = {0.0};
    double solution[64];
    c_ode_settings settings;
    c_ode_event_solver_settings esettings;
    c_ode_event_statistics estats;
    c_ode_results results = NULL;
    int rows, columns, i;

    c_default_ode_settings(&settings);
    c_default_ode_event_solver_settings(&esettings);
    settings.absolute_tolerance = 1.0e-10;
    settings.relative_tolerance = 1.0e-10;
    settings.maximum_step_size = 0.1;
    event_callback_calls = 0;
    callback_errors = 0;

    check(solve(2, 1, 2, unit_rate, x, initial, &settings, &esettings,
        two_events, &results, &estats) == C_NO_ERROR,
        "event-enabled integration succeeds");
    check(results != NULL, "event-enabled integration returns results");
    check(callback_errors == 0, "event callback receives declared dimensions");
    check(event_callback_calls > 0, "event callback is evaluated");
    check(estats.event_index == 1, "statistics report zero-based earliest event");
    check(estats.function_evaluations > 0, "event evaluation count is positive");
    check(estats.iteration_count > 0, "event localization takes iterations");
    check(c_get_ode_results_size(results, &rows, &columns) == C_NO_ERROR,
        "get event result dimensions");
    check(rows >= 2 && rows <= 64 && columns == 2,
        "event result has time and state columns");
    check(c_get_ode_results(results, rows, columns, solution, rows) == C_NO_ERROR,
        "copy event result");
    for (i = 0; i < rows; ++i) {
        check(solution[i] <= 0.3 + 1.0e-8,
            "event result does not advance past earliest event");
    }
    check_close(solution[rows - 1], 0.3, 1.0e-8,
        "last result time is the earliest event");
    check_close(solution[rows + rows - 1], 0.3, 1.0e-8,
        "last result state is evaluated at the event");
    c_free_ode_results(results);
}

static void test_no_event_fallback(event_integrator_function solve, int nevents,
    c_ode_event_function efcn)
{
    const double x[] = {0.0, 1.0};
    const double initial[] = {0.0};
    double solution[64];
    c_ode_settings settings;
    c_ode_event_solver_settings esettings;
    c_ode_event_statistics estats;
    c_ode_results results = NULL;
    int rows, columns;

    c_default_ode_settings(&settings);
    c_default_ode_event_solver_settings(&esettings);
    settings.absolute_tolerance = 1.0e-10;
    settings.relative_tolerance = 1.0e-10;
    event_callback_calls = 0;

    callback_errors = 0;
    check(solve(2, 1, nevents, unit_rate, x, initial, &settings, &esettings,
        efcn, &results, &estats) == C_NO_ERROR,
        "disabled-event fallback integration succeeds");
    check(results != NULL, "fallback integration returns results");
    check(event_callback_calls == 0, "disabled event callback is not invoked");
    check(callback_errors == 0, "disabled event callback is not invoked");
    check(estats.event_index == -1, "no event uses the sentinel event index");
    check(c_get_ode_results_size(results, &rows, &columns) == C_NO_ERROR,
        "get fallback result dimensions");
    check(rows >= 2 && rows <= 64 && columns == 2,
        "fallback returns the requested interval");
    check(c_get_ode_results(results, rows, columns, solution, rows) == C_NO_ERROR,
        "copy fallback result");
    check_close(solution[rows - 1], 1.0, 1.0e-12,
        "fallback reaches final requested time");
    check_close(solution[2 * rows - 1], 1.0, 1.0e-8,
        "fallback reaches final state");
    c_free_ode_results(results);
}

int main(int argc, char **argv)
{
    size_t i;
    int matched = 0;

    test_context = "arguments";
    check(argc <= 2, "expected at most one test name");
    for (i = 0; i < sizeof(integrators) / sizeof(integrators[0]); ++i) {
        if (argc == 2 && strcmp(argv[1], integrators[i].name) != 0) continue;
        matched = 1;
        test_context = integrators[i].name;
        test_event_detection(integrators[i].solve);
        test_no_event_fallback(integrators[i].solve, 2, NULL);
        test_no_event_fallback(integrators[i].solve, 0, two_events);
        test_no_event_fallback(integrators[i].solve, -1, two_events);
    }
    check(matched, "unknown test name");
    return EXIT_SUCCESS;
}
