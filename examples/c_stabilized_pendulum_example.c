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

const double body_mass = 2.0;
const double length = 1.5;

void force(int n, double t, const double *x, double *frc);
void physical_mass(int n, double t, const double *x, double *m, int ldm);
void constraints(int nc, int nstates, double t, const double *x, double *g);
void gradient(int m, int n, double t, const double *x, double *grad, int ldg);
void curvature(int nc, int nstates, double t, const double *q, 
    const double *direction, const double *v, double *crv);

int main()
{
    // Local Variables
    c_ggl_container model;
    double ic[] = {length, 0.0, 0.0, 0.0, 0.0, 0.0};
    double t[] = {0.0, 10.0};
    c_ode_results rst;
    c_ode_settings settings;
    int flag;

    /* NOTE:
    * The state vector for this example is a 6 component vector:
    * 1. x motion
    * 2. y motion
    * 3. x velocity
    * 4. y velocity
    * 5. lambda
    * 6. mu
    */

    // Set up the model
    model.alpha = 2.0;
    model.beta = 3.0;
    model.constraint_count = 1;
    model.constraints = constraints;
    model.curvature = curvature;
    model.force = force;
    model.gradient = gradient;
    model.mass = physical_mass;
    model.jacobian = NULL;

    // Set up the integrator
    c_default_ode_settings(&settings);

    // Solve
    flag = c_dae_bdf(&model, 2, 6, t, ic, &settings, &rst);
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Constrained integration failed: %d\n", flag);
        return flag;
    }

    // Write the results to file
    flag = c_write_ode_results("pendulum_example.txt", rst, ',');
    if (flag != C_NO_ERROR) {
        fprintf(stderr, "Writing results failed: %d\n", flag);
    }

    // Clean up
    c_free_ode_results(rst);
    rst = NULL;

    // End
    return flag;
}


void force(int n, double t, const double *x, double *frc)
{
    frc[0] = 0.0;
    frc[1] = -body_mass * 9.81;
}

void physical_mass(int n, double t, const double *x, double *m, int ldm)
{
    m[0] = body_mass;
    m[1] = 0.0;
    m[2] = 0.0;
    m[3] = body_mass;
}

void constraints(int nc, int nstates, double t, const double *x, double *g)
{
    g[0] = x[0] * x[0] + x[1] * x[1] - length * length;
}

void gradient(int m, int n, double t, const double *x, double *grad, int ldg)
{
    int i;
    for (i = 0; i < n; ++i) grad[i] = 2.0 * x[i];
}

void curvature(int nc, int nstates, double t, const double *q, 
    const double *direction, const double *v, double *crv)
{
    crv[0] = 2.0 * (direction[0] * v[0] + direction[1] * v[1]);
}
