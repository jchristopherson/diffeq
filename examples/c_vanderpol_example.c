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
#include <stdlib.h>
#include "diffeq.h"

void vanderpol(int n, double x, const double *y, double *dydx);

int main()
{
    // Local Variables
    const int npts = 2;
    const int nvars = 2;
    double x[] = {0.0, 50.0};
    double ic[] = {2.0, 0.0};
    int m, n, flag;
    double *sol;
    c_ode_results rst;
    c_ode_settings settings;

    // Use the default integrator settings
    c_default_ode_settings(&settings);

    // Integrate using a 5th order Tsitouras integrator.
    flag = c_tsitouras(npts, nvars, vanderpol, x, ic, &settings, &rst);
    if (flag != C_NO_ERROR) {
        printf("There was an error encountered during integration.  Flag: %i\n", flag);
        return flag;
    }

    // Utilize the c_ode_results object to extract all the solution points as
    // a matrix where the first column contains the values of x at which the
    // solution was computed.  The remaining 2 columns contain the solution
    // points for each state variable.  Utilize c_get_ode_results_size to
    // get the size of the matrix required to hold the results.
    flag = c_get_ode_results_size(rst, &m, &n);
    if (flag != C_NO_ERROR)
    {
        printf("There was an error retrieving the size of the results. Flag: %i\n", flag);
        return flag;
    }
    sol = (double*)malloc((size_t)(m * n * sizeof(double)));
    if (!sol)
    {
        printf("There was an issue allocating memory for the results matrix.\n");
        c_free_ode_results(rst);
        return -1;
    }

    // Now that we have our matrix allocated, we can retrieve the results.
    flag = c_get_ode_results(rst, m, n, sol, m);
    if (flag != C_NO_ERROR)
    {
        printf("There was an issue copying the results matrix.  Flag: %i\n", flag);
        free(sol);
        c_free_ode_results(rst);
        return flag;
    }

    // An alternative to copying the results to a matrix is to write the results
    // directly to a delimited text file.  This can be accomplished as follows.
    flag = c_write_ode_results("example file.txt", rst, ',');
    if (flag != C_NO_ERROR)
    {
        printf("There was an issue writing the results to file.  Flag: %i\n", flag);
        free(sol);
        c_free_ode_results(rst);
        return flag;
    }

    // Clean up
    free(sol);
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
