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

#include <stdio.h>
#include "diffeq.h"

void vanderpol(int n, double x, const double *y, double *dydx);

int main()
{
    // Local Variables
    double x[] = {0.0, 50.0};
    double ic[] = {2.0, 0.0};
    int n, npts;
    c_ode_results rst;
    c_ode_results_statistics stats;
    c_ode_settings settings;

    // Use default settings
    c_default_ode_settings(&settings);

    // RK23
    c_runge_kutta_23(2, 2, vanderpol, x, ic, &settings, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    
    printf("RUNGE_KUTTA_23:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // RK45
    c_runge_kutta_45(2, 2, vanderpol, x, ic, &settings, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    
    printf("RUNGE_KUTTA_45:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // RK853
    c_runge_kutta_853(2, 2, vanderpol, x, ic, &settings, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    
    printf("RUNGE_KUTTA_853:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // ROSENBROCK
    c_rosenbrock(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("ROSENBROCK:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Jacobian Evaluations\n", stats.jacobian_evaluations);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    printf("\t%i Matrix Factorizations\n", stats.matrix_factorizations);
    c_free_ode_results(rst);

    // RK45 W/ PI
    settings.beta = 0.1;
    c_runge_kutta_45(2, 2, vanderpol, x, ic, &settings, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("RUNGE_KUTTA_45 w/ PI Controller:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // KC4
    c_kennedy_carpenter_4(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("KC4:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Jacobian Evaluations\n", stats.jacobian_evaluations);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    printf("\t%i Matrix Factorizations\n", stats.matrix_factorizations);
    c_free_ode_results(rst);

    // KC5
    c_kennedy_carpenter_5(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("KC5:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Jacobian Evaluations\n", stats.jacobian_evaluations);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    printf("\t%i Matrix Factorizations\n", stats.matrix_factorizations);
    c_free_ode_results(rst);

    // TSITOURAS
    c_tsitouras(2, 2, vanderpol, x, ic, &settings, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    
    printf("TSITOURAS 4/5:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // BDF
    c_bdf(2, 2, vanderpol, x, ic, &settings, NULL, NULL, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("BDF:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Jacobian Evaluations\n", stats.jacobian_evaluations);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    printf("\t%i Matrix Factorizations\n", stats.matrix_factorizations);
    c_free_ode_results(rst);

    // ADAMS
    c_adams(2, 2, vanderpol, x, ic, &settings, NULL, &rst);
    c_get_ode_results_statistics(rst, &stats);
    c_get_ode_results_size(rst, &npts, &n);
    settings.beta = 0.0;    // reset to avoid using PI with other solvers
    
    printf("ADAMS:\n");
    printf("\t%i Solution Points\n", npts);
    printf("\t%i Successful Steps\n", stats.successful_steps);
    printf("\t%i Failed Steps\n", stats.failed_steps);
    printf("\t%i Function Evaluations\n", stats.function_evaluations);
    c_free_ode_results(rst);

    // End
    return 0;
}


void vanderpol(int n, double x, const double *y, double *dydx)
{
    const double mu = 5.0;
    dydx[0] = y[1];
    dydx[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
}
