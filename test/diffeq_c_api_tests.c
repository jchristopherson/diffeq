#include "diffeq.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { expected_rows = 5, expected_cols = 3, leading_dimension = 7 };
static const double sample_times[] = {0.0, 0.25, 0.5, 0.75, 1.0};
static const double expected_exp[] = {
    1.0, 1.2840254167, 1.6487212707, 2.1170000166, 2.7182818285
};
static const double expected_linear[] = {-1.0, -0.5, 0.0, 0.5, 1.0};

static void check(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "C API test failed: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void exponential_system(int n, double x, const double *y, double *dydx)
{
    (void)x;
    check(n == 2, "callback receives expected state count");
    dydx[0] = y[0];
    dydx[1] = 2.0;
}

static void check_written_row(char *line, const double expected[3])
{
    char *cursor = line;
    int column;

    for (column = 0; column < 3; ++column) {
        char *end;
        double value = strtod(cursor, &end);
        check(end != cursor, "parse written result value");
        check(fabs(value - expected[column]) < 1.0e-5,
            "written result value");
        cursor = end;
        if (column < 2) {
            check(*cursor == ';', "written result delimiter");
            ++cursor;
        }
    }
    check(*cursor == '\n' || *cursor == '\0', "written result row terminator");
}

static void test_default_settings(void)
{
    c_ode_settings settings;

    c_default_ode_settings(&settings);
    check(settings.mass_matrix_constant, "default constant mass matrix");
    check(settings.minimum_step_size == 10.0 * DBL_EPSILON,
        "default minimum step size");
    check(settings.maximum_step_size == DBL_MAX, "default maximum step size");
    check(settings.step_limit == 1000000, "default step limit");
    check(settings.allow_overshoot, "default overshoot setting");
    check(settings.beta == 0.0, "default step size controller");
    check(settings.absolute_tolerance == 1.0e-6, "default absolute tolerance");
    check(settings.relative_tolerance == 1.0e-6, "default relative tolerance");
}

static void test_null_results(void)
{
    c_ode_results_statistics stats;
    int rows = -1;
    int cols = -1;
    double matrix[1];

    c_free_ode_results(NULL);
    check(c_get_ode_results_statistics(NULL, &stats) == C_NULL_POINTER_ERROR,
        "statistics rejects null result");
    check(c_get_ode_results_size(NULL, &rows, &cols) == C_NULL_POINTER_ERROR,
        "size rejects null result");
    check(rows == 0 && cols == 0, "null size initializes dimensions");
    check(c_get_ode_results(NULL, 0, 0, matrix, 0) == C_NULL_POINTER_ERROR,
        "matrix getter rejects null result");
    check(c_write_ode_results("diffeq_c_api_null.txt", NULL, ';') ==
        C_NULL_POINTER_ERROR, "writer rejects null result");
}

static c_ode_results create_test_results(void)
{
    const double initial[] = {1.0, -1.0};
    c_ode_settings settings;
    c_ode_results results = NULL;

    c_default_ode_settings(&settings);
    check(c_runge_kutta_23(expected_rows, 2, exponential_system, sample_times,
        initial, &settings, &results) == C_NO_ERROR, "create result fixture");
    check(results != NULL, "integration returns result handle");
    return results;
}

static void test_result_matrix(c_ode_results results)
{
    double matrix[leading_dimension * expected_cols];
    c_ode_results_statistics stats;
    int rows, cols, index;

    check(c_get_ode_results_size(results, &rows, &cols) == C_NO_ERROR,
        "get result dimensions");
    check(rows == expected_rows && cols == expected_cols, "result dimensions");

    check(c_get_ode_results_statistics(results, &stats) == C_NO_ERROR,
        "get integration statistics");
    check(stats.successful_steps > 0, "successful step count");
    check(stats.failed_steps >= 0, "failed step count");
    check(stats.function_evaluations > 0, "function evaluation count");

    check(c_get_ode_results(results, rows - 1, cols, matrix,
        leading_dimension) == C_MATRIX_ROW_SIZE_ERROR,
        "matrix getter validates row count");
    check(c_get_ode_results(results, rows, cols - 1, matrix,
        leading_dimension) == C_MATRIX_COLUMN_SIZE_ERROR,
        "matrix getter validates column count");
    check(c_get_ode_results(results, rows, cols, matrix, rows - 1) ==
        C_MATRIX_STRIDE_ERROR, "matrix getter validates leading dimension");

    for (index = 0; index < leading_dimension * expected_cols; ++index) {
        matrix[index] = -12345.0;
    }
    check(c_get_ode_results(results, rows, cols, matrix, leading_dimension) ==
        C_NO_ERROR, "copy result matrix");
    for (index = 0; index < expected_rows; ++index) {
        check(fabs(matrix[index] - sample_times[index]) < 1.0e-12,
            "independent variable column");
        check(fabs(matrix[index + leading_dimension] - expected_exp[index]) <
            1.0e-5, "exponential state column");
        check(fabs(matrix[index + 2 * leading_dimension] -
            expected_linear[index]) < 1.0e-12, "linear state column");
    }
    for (index = expected_rows; index < leading_dimension; ++index) {
        check(matrix[index] == -12345.0, "first column stride padding");
        check(matrix[index + leading_dimension] == -12345.0,
            "second column stride padding");
        check(matrix[index + 2 * leading_dimension] == -12345.0,
            "third column stride padding");
    }
}

static void test_result_writer(c_ode_results results)
{
    const char *filename = "diffeq_c_api_results.txt";
    FILE *file;
    char line[256];
    int index;

    check(c_write_ode_results(filename, results, ';') == C_NO_ERROR,
        "write integration results");
    file = fopen(filename, "r");
    check(file != NULL, "open written results");
    for (index = 0; index < expected_rows; ++index) {
        check(fgets(line, sizeof(line), file) != NULL, "read written row");
        {
            double expected[3] = {sample_times[index], expected_exp[index],
                expected_linear[index]};
            check_written_row(line, expected);
        }
    }
    check(fgets(line, sizeof(line), file) == NULL, "writer emits only results");
    check(fclose(file) == 0, "close written results");
    check(remove(filename) == 0, "remove written results");
}

int main(void)
{
    c_ode_results results;

    test_default_settings();
    test_null_results();
    results = create_test_results();
    test_result_matrix(results);
    test_result_writer(results);
    c_free_ode_results(results);
    return EXIT_SUCCESS;
}
