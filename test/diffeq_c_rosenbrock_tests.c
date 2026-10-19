#include "diffeq.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int use_jacobian;
static int use_mass;
static int jacobian_calls;
static int mass_calls;

static void check(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "C Rosenbrock test failed [Jacobian=%d, mass=%d]: %s\n",
            use_jacobian, use_mass, message);
        exit(EXIT_FAILURE);
    }
}

static void coupled_system(int n, double x, const double *y, double *dydx)
{
    (void)x;
    check(n == 2, "coupled RHS state count");
    dydx[0] = -y[0] + y[1];
    dydx[1] = -2.0 * y[1];
}

static void coupled_jacobian(int m, int n, double x, const double *y,
    double *jac, int ldj)
{
    (void)x;
    (void)y;
    check(m == 2 && n == 2 && ldj >= m, "Jacobian dimensions and stride");
    ++jacobian_calls;
    jac[0] = -1.0;
    jac[1] = 0.0;
    jac[ldj] = 1.0;
    jac[ldj + 1] = -2.0;
}

/* M = [2 1; 0 3], so M*y' = [-2*y0 - 2*y1; -6*y1]. */
static void mass_system(int n, double x, const double *y, double *dydx)
{
    (void)x;
    check(n == 2, "mass RHS state count");
    dydx[0] = -2.0 * y[0] - 2.0 * y[1];
    dydx[1] = -6.0 * y[1];
}

static void mass_jacobian(int m, int n, double x, const double *y,
    double *jac, int ldj)
{
    (void)x;
    (void)y;
    check(m == 2 && n == 2 && ldj >= m, "mass Jacobian dimensions and stride");
    ++jacobian_calls;
    jac[0] = -2.0;
    jac[1] = 0.0;
    jac[ldj] = -2.0;
    jac[ldj + 1] = -6.0;
}

static void mass_matrix(int n, double x, const double *y, double *mass, int ldm)
{
    (void)x;
    (void)y;
    check(n == 2 && ldm >= n, "mass matrix dimensions and stride");
    ++mass_calls;
    mass[0] = 2.0;
    mass[1] = 0.0;
    mass[ldm] = 1.0;
    mass[ldm + 1] = 3.0;
}

static void test_optional_callbacks(void)
{
    const double x[] = {0.0, 0.1, 0.4, 1.0};
    const double initial[] = {1.0, 1.0};
    c_ode_settings settings;
    c_ode_results results = NULL;
    c_ode_results_statistics stats;
    c_ode_jacobian jac = use_mass ? mass_jacobian : coupled_jacobian;
    double matrix[12];
    int rows, cols, i;

    c_default_ode_settings(&settings);
    settings.absolute_tolerance = 1.0e-9;
    settings.relative_tolerance = 1.0e-9;
    jacobian_calls = mass_calls = 0;
    check(c_rosenbrock(4, 2, use_mass ? mass_system : coupled_system,
        x, initial, &settings, use_jacobian ? jac : NULL,
        use_mass ? mass_matrix : NULL, &results) == C_NO_ERROR,
        "solve callback system");
    check(results != NULL, "callback system returns a result handle");
    check(c_get_ode_results_size(results, &rows, &cols) == C_NO_ERROR,
        "callback result size");
    check(rows == 4 && cols == 3, "callback result dimensions");
    check(c_get_ode_results(results, rows, cols, matrix, rows) == C_NO_ERROR,
        "copy callback solution");
    for (i = 0; i < rows; ++i) {
        double e1 = exp(-x[i]);
        double e2 = exp(-2.0 * x[i]);
        check(fabs(matrix[i] - x[i]) < 1.0e-12, "callback solution time");
        check(fabs(matrix[i + rows] - (use_mass ? e1 : 2.0 * e1 - e2)) <
            1.0e-6, "coupled first state");
        check(fabs(matrix[i + 2 * rows] - e2) < 1.0e-6,
            "coupled second state");
    }
    check(c_get_ode_results_statistics(results, &stats) == C_NO_ERROR,
        "callback statistics");
    check(stats.jacobian_evaluations > 0 && stats.matrix_factorizations > 0,
        "Rosenbrock performs Jacobian evaluations and factorizations");
    check(use_jacobian ? jacobian_calls > 0 : jacobian_calls == 0,
        "optional Jacobian callback dispatch");
    check(use_mass ? mass_calls > 0 : mass_calls == 0,
        "optional mass matrix callback dispatch");
    c_free_ode_results(results);
}

int main(void)
{
    for (use_jacobian = 0; use_jacobian <= 1; ++use_jacobian) {
        for (use_mass = 0; use_mass <= 1; ++use_mass) {
            test_optional_callbacks();
        }
    }
    return EXIT_SUCCESS;
}
