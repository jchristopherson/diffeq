#include "diffeq.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int (*integrator_function)(int, int, c_ode, const double *,
    const double *, const c_ode_settings *, c_ode_results *);

static const char *test_context;
static int rhs_calls;
static int reverse_case;
static int overshoot_case;

static void check(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "C integrator test failed [%s, reverse=%d, overshoot=%d]: %s\n",
            test_context, reverse_case, overshoot_case, message);
        exit(EXIT_FAILURE);
    }
}

static void check_close(double actual, double expected, const char *message)
{
    if (!(fabs(actual - expected) < 1.0e-6)) {
        fprintf(stderr, "Expected %.17g, received %.17g\n", expected, actual);
        check(0, message);
    }
}

static int rosenbrock_default(int nx, int nstates, c_ode fcn, const double *x,
    const double *initial, const c_ode_settings *settings, c_ode_results *output)
{
    return c_rosenbrock(nx, nstates, fcn, x, initial, settings, NULL, NULL,
        output);
}

static const struct {
    const char *name;
    integrator_function solve;
} integrators[] = {
    {"rk23", c_runge_kutta_23},
    {"rk45", c_runge_kutta_45},
    {"rk853", c_runge_kutta_853},
    {"tsitouras", c_tsitouras},
    {"adams", c_adams},
    {"rosenbrock", rosenbrock_default}
};

static void exponential_system(int n, double x, const double *y, double *dydx)
{
    check(n == 2, "RHS state count");
    check(isfinite(x), "RHS independent variable");
    ++rhs_calls;
    dydx[0] = y[0];
    dydx[1] = 2.0 * x;
}

static void test_invalid_input(integrator_function solve)
{
    const double x[] = {0.0, 1.0};
    const double initial[] = {1.0, -1.0};
    c_ode_settings settings;
    c_ode_results results = NULL;

    c_default_ode_settings(&settings);
    rhs_calls = 0;
    check(solve(2, 2, NULL, x, initial, &settings, &results) ==
        C_NULL_POINTER_ERROR, "null RHS callback");
    check(solve(1, 2, exponential_system, x, initial, &settings, &results) ==
        C_INVALID_INPUT_ERROR, "short grid");
    check(solve(0, 2, exponential_system, x, initial, &settings, &results) ==
        C_INVALID_INPUT_ERROR, "empty grid");
    check(solve(2, 0, exponential_system, x, initial, &settings, &results) ==
        C_INVALID_INPUT_ERROR, "zero states");
    check(solve(2, -1, exponential_system, x, initial, &settings, &results) ==
        C_INVALID_INPUT_ERROR, "negative state count");
    check(rhs_calls == 0, "invalid input does not invoke RHS");
}

static void test_solution(integrator_function solve, int reverse,
    int allow_overshoot)
{
    double x[5];
    double initial[2];
    double matrix[15];
    c_ode_settings settings;
    c_ode_results results = NULL;
    c_ode_results_statistics stats;
    int rows, cols, i;

    reverse_case = reverse;
    overshoot_case = allow_overshoot;
    for (i = 0; i < 5; ++i) {
        x[i] = reverse ? 1.0 - 0.25 * i : 0.25 * i;
    }
    initial[0] = exp(x[0]);
    initial[1] = x[0] * x[0] - 1.0;
    c_default_ode_settings(&settings);
    settings.allow_overshoot = allow_overshoot != 0;
    settings.absolute_tolerance = 1.0e-9;
    settings.relative_tolerance = 1.0e-9;
    settings.maximum_step_size = 0.1;
    rhs_calls = 0;

    check(solve(5, 2, exponential_system, x, initial, &settings, &results) ==
        C_NO_ERROR, "solve analytical system");
    check(results != NULL, "non-null result handle");
    check(c_get_ode_results_size(results, &rows, &cols) == C_NO_ERROR,
        "get result size");
    check(rows == 5 && cols == 3, "requested-grid result dimensions");
    check(c_get_ode_results(results, rows, cols, matrix, rows) == C_NO_ERROR,
        "copy solution");
    for (i = 0; i < rows; ++i) {
        check(fabs(matrix[i] - x[i]) < 1.0e-12, "requested solution time");
        check_close(matrix[i + rows], exp(x[i]),
            "state-dependent exponential solution");
        check_close(matrix[i + 2 * rows], x[i] * x[i] - 1.0,
            "time-dependent polynomial solution");
    }
    check(c_get_ode_results_statistics(results, &stats) == C_NO_ERROR,
        "get statistics");
    check(stats.successful_steps > 0 && stats.failed_steps >= 0,
        "step statistics");
    check(rhs_calls > 0 && stats.function_evaluations == rhs_calls,
        "RHS callback count agrees with statistics");
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
        test_invalid_input(integrators[i].solve);
        test_solution(integrators[i].solve, 0, 1);
        test_solution(integrators[i].solve, 0, 0);
        test_solution(integrators[i].solve, 1, 1);
        test_solution(integrators[i].solve, 1, 0);
    }
    check(matched, "unknown test name");
    return EXIT_SUCCESS;
}
