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
    const double *, double *, const c_ode_settings *, c_ode_results *);

static const char *test_name;
static int force_calls;
static int mass_calls;
static int constraint_calls;
static int gradient_calls;
static int curvature_calls;

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
    c_ode_results *output)
{
    return c_dae_rosenbrock(model, nx, nstates, x, initial, settings, output);
}

static int dae_kc4(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output)
{
    return c_dae_kennedy_carpenter_4(model, nx, nstates, x, initial, settings,
        output);
}

static int dae_kc5(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output)
{
    return c_dae_kennedy_carpenter_5(model, nx, nstates, x, initial, settings,
        output);
}

static int dae_bdf(int nx, int nstates, const c_ggl_container *model,
    const double *x, double *initial, const c_ode_settings *settings,
    c_ode_results *output)
{
    return c_dae_bdf(model, nx, nstates, x, initial, settings, output);
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

static void test_solution(dae_integrator solve)
{
    const double times[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    double initial[state_count] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double matrix[5 * (state_count + 1)];
    c_ggl_container model = {0};
    c_ode_settings settings;
    c_ode_results results = NULL;
    c_ode_results_statistics stats;
    int rows, cols, i;

    model.force = oscillator_force;
    model.mass = identity_mass;
    model.constraints = constrain_y;
    model.gradient = constrain_y_gradient;
    model.curvature = linear_constraint_curvature;
    model.alpha = 2.0;
    model.beta = 3.0;
    model.constraint_count = constraint_count;

    c_default_ode_settings(&settings);
    settings.absolute_tolerance = 1.0e-9;
    settings.relative_tolerance = 1.0e-9;
    settings.maximum_step_size = 0.05;
    force_calls = mass_calls = constraint_calls = gradient_calls =
        curvature_calls = 0;

    check(solve(5, state_count, &model, times, initial, &settings, &results) ==
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
    c_free_ode_results(results);
}

int main(int argc, char **argv)
{
    size_t i;

    check(argc == 2, "expected one DAE integrator name");
    for (i = 0; i < sizeof(integrators) / sizeof(integrators[0]); ++i) {
        if (strcmp(argv[1], integrators[i].name) == 0) {
            test_name = integrators[i].name;
            test_solution(integrators[i].solve);
            return EXIT_SUCCESS;
        }
    }
    test_name = "arguments";
    check(0, "unknown DAE integrator name");
    return EXIT_FAILURE;
}
