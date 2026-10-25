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

#define INDEX(i, j, ld) ((j) * (ld) + (i))

void ball(int n, double t, const double *y, double *dydt);
void ball_ground_constraint(int neqn, int nstats, double t, const double *y, double *f);
void clean_up(c_ode_results r1, c_ode_results r2, double *s1, double *s2);
void print_event_statistics(int event, const c_ode_event_statistics *x);
void print_matrix(int m, int n, const double *x, int ldx);

int main()
{
    // Local Variables
    const int nt = 2;
    const int neqn = 2;
    const int nevents = 1;
    const double k = 0.7;
    int flag, m1, m2, n;
    double t[] = {0.0, 10.0};
    double ic[] = {5.0, 0.0};
    c_ode_settings settings;
    c_ode_event_solver_settings esettings;
    c_ode_results r1, r2;
    c_ode_event_statistics e1, e2;
    double *s1, *s2;

    // Get the default settings
    c_default_ode_settings(&settings);
    c_default_ode_event_solver_settings(&esettings);

    // Solve for the first bounce
    flag = c_tsitouras_events(nt, neqn, nevents, ball, t, ic, &settings, 
        &esettings, ball_ground_constraint, &r1, &e1);
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Integration failed: %d\n", flag);
        return flag;
    }
    
    c_get_ode_results_size(r1, &m1, &n);
    s1 = (double*)malloc((size_t)(m1 * n * sizeof(double)));
    if (!s1) {
        fprintf(stderr, "Memory allocation failed:\n");
        clean_up(r1, r2, s1, s2);
        return EXIT_FAILURE;
    }
    flag = c_get_ode_results(r1, m1, n, s1, m1);
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Error retrieving results from set 1: %d\n", flag);
        clean_up(r1, r2, s1, s2);
        return flag;
    }
    print_event_statistics(1, &e1);

    flag = c_write_ode_results("Event Bounce 1.txt", r1, ',');
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Error writing results file from set 1: %d\n", flag);
        clean_up(r1, r2, s1, s2);
        return flag;
    }

    // Now, handle the second bounce.  Reset the initial conditions and 
    // starting time point
    t[0] = s1[INDEX(m1 - 1, 0, m1)];
    ic[0] = 0.0;
    ic[1] = -k * s1[INDEX(m1 - 1, 2, m1)];
    flag = c_tsitouras_events(nt, neqn, nevents, ball, t, ic, &settings,
        &esettings, ball_ground_constraint, &r2, &e2);

    c_get_ode_results_size(r2, &m2, &n);
    s2 = (double*)malloc((size_t)(m2 * n * sizeof(double)));
    if (!s2) {
        fprintf(stderr, "Memory allocation failed:\n");
        clean_up(r1, r2, s1, s2);
        return EXIT_FAILURE;
    }
    flag = c_get_ode_results(r2, m2, n, s2, m2);
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Error retrieving results from set 2: %d\n", flag);
        clean_up(r1, r2, s1, s2);
        return flag;
    }
    print_event_statistics(2, &e2);

    flag = c_write_ode_results("Event Bounce 2.txt", r2, ',');
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Error writing results file from set 2: %d\n", flag);
        clean_up(r1, r2, s1, s2);
        return flag;
    }

    // Print out the results
    printf("Time [s]\tHeight [m]\tVelocity [m/s]\n");
    print_matrix(m1, n, s1, m1);

    printf("\n");
    print_matrix(m2, n, s2, m2);

    // End
    clean_up(r1, r2, s1, s2);
    return 0;
}


void ball(int n, double t, const double *y, double *dydt) 
{
    const double g = 9.81;
    dydt[0] = y[1];
    dydt[1] = -g;
}

void ball_ground_constraint(int neqn, int nstats, double t, const double *y, double *f)
{
    f[0] = y[0];
}

void clean_up(c_ode_results r1, c_ode_results r2, double *s1, double *s2)
{
    if (r1) c_free_ode_results(r1);
    if (r2) c_free_ode_results(r2);
    if (s1) free(s1);
    if (s2) free(s2);
}

void print_event_statistics(int event, const c_ode_event_statistics *x)
{
    printf("Event %d Statistics:\n", event);
    printf("\tEvent Index: %d\n", x->event_index);
    printf("\tEvent Function Evaluations: %d\n", x->function_evaluations);
    printf("\tEvent Solver Iterations: %d \n", x->iteration_count);
}

void print_matrix(int m, int n, const double *x, int ldx)
{
    int i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            printf("%9.6f\t", x[INDEX(i, j, ldx)]);
        }
        printf("\n");
    }
}