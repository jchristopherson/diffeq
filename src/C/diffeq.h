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

#ifndef DIFFEQ_H_
#define DIFFEQ_H_

#include <stdbool.h>

#define C_NO_ERROR 0
#define C_NULL_POINTER_ERROR 1
#define C_MATRIX_ROW_SIZE_ERROR 2
#define C_MATRIX_COLUMN_SIZE_ERROR 3
#define C_MATRIX_STRIDE_ERROR 4
#define C_INVALID_INPUT_ERROR 5

/**
 * @brief ODE callback routine.
 * 
 * @param n The number of equations being integrated.
 * @param x The current value of the independent variable.
 * @param y The N-element state vector.
 * @param dydx An N-element vector where the state vector derivatives are written.
 */
typedef void (*c_ode)(int n, double x, const double *y, double *dydx);

/**
 * @brief ODE Jacobian callback routine.
 * 
 * @param m The number of equations.
 * @param n The state vector length.
 * @param x The current value of the independent variable.
 * @param y The N-element state vector.
 * @param jac The M-by-N Jacobian matrix stored in column-major format.
 * @param ldj The leading dimension of @p jac.  This value must be at least @p m.
 */
typedef void (*c_ode_jacobian)(int m, int n, double x, const double *y, double *jac, int ldj);

/**
 * @brief ODE mass matrix callback routine.
 * 
 * @param n The state vector length.
 * @param x The current value of the independent variable.
 * @param y The N-element state vector.
 * @param m The N-by-N mass matrix stored in column-major format.
 * @param ldm The leading dimension of @p m.  This value must be at least @p n.
 */
typedef void (*c_ode_mass_matrix)(int n, double x, const double *y, double *m, int ldm);

/**
 * @brief Evaluates algebraic constraint residuals \(g(x,y)\).
 *
 * For \(N\) input state components and \(M\) constraints,
 * the callback represents a mapping
 * \[
 * g : \mathbb{R}\times\mathbb{R}^{N}
 * \longrightarrow \mathbb{R}^{M}.
 * \]
 * A consistent state satisfies \(g(x,y)=0\). Return the
 * residuals themselves, not their derivatives, a Boolean result,
 * or a projected state. A nonzero residual is a valid evaluation
 * at an inconsistent state, rather than a callback failure.
 *
 * The caller supplies the input and output array sizes. The
 * implementation must interpret the same state ordering as its
 * model and write every output element. This interface performs
 * no dimension, rank, or initial-consistency checks. Residual
 * ordering and physical units are defined by the model; any
 * subsequent consistency tolerance must respect that scaling.
 * 
 * @param nconstraints The number of constraint equations.
 * @param nstates The number of state variables.
 * @param x The current value of the independent variable.
 * @param y The current state vector.
 * @param g The constraint vector.
 */
typedef void (*c_dae_algebraic_constraints)(int nconstraints, int nstates,
    double x, const double *y, double *g);

/**
 * @brief Evaluates \(DG(q)[\mathrm{direction}]v\), without including
 * \(G\dot v\). Each component is
 * \(\mathrm{direction}^T\nabla^2 g_i(q)v\).
 * Assign every output element. For linear constraints return
 * zeros. The callback must support off-manifold coordinates and
 * arbitrary directions used during nonlinear/Jacobian evaluation.
 * 
 * @param nconstraints The number of constraint equations.
 * @param nstates The number of state variables.
 * @param x The current value of the independent variable.
 * @param q The current coordinate vector.
 * @param direction The coordinate derivative \(\dot q\), including the 
 * kinematic multiplier correction.
 * @param v The physical velocity components.
 * @param curvature The curvature terms in the order of teh constraint callback.
 */
typedef void (*c_dae_constraint_curvature)(int nconstraints, int nstates, 
    double x, const double *q, const double *direction, const double *v,
    double *curvature);

/**
 * @brief A container of performance statistics for an integrator.
 */
typedef struct {
    /** The number of successful steps taken by the integrator. */
    int successful_steps;
    /** The number of unsuccessful steps taken by the integrator. */
    int failed_steps;
    /** The number of Jacobian matrix evaluations. */
    int jacobian_evaluations;
    /** The number of function evaluations. */
    int function_evaluations;
    /** The number of matrix factorizations. */
    int matrix_factorizations;
} c_ode_results_statistics;

/** @brief A container for storing results from the integrator. Use c_free_ode_results
 * to release resources once finished using this object.
*/
typedef void* c_ode_results;

/**
 * @brief A container for ODE solver settings.
 */
typedef struct {
    /** 
     * True if the mass matrix is constant regardless of state, or false if
     * the mass matrix is state dependent.
     */
    bool mass_matrix_constant;
    /** The minimum allowed integrator step size. */
    double minimum_step_size;
    /** The maximum allowed integrator step size. */
    double maximum_step_size;
    /** The limit on the maximum total number of integration steps. */
    int step_limit;
    /**
     * True if the solver is allowed to overshoot the final value in the 
     * integration range and interpolate back; else, false if the solver must
     * terminate on the final value.
     */
    bool allow_overshoot;
    /** 
     * PI step size controller exponent.  Set to zero to turn off PI step size
     * control.
     */
    double beta;
    /** The absolute tolerance applied to each equation. */
    double absolute_tolerance;
    /** The relative tolerance applied to each equation. */
    double relative_tolerance;
    /** The maximum number of Newton iterations allowed per implicit solve
     * (BDF and Kennedy-Carpenter). */
    int newton_step_limit;
    /** The convergence tolerance for Newton iteration (BDF and
     * Kennedy-Carpenter). */
    double newton_tolerance;
} c_ode_settings;

/**
 * @brief Describes a holonomically constrained mechanical model to be solved
 * with the stabilized index-1 Gear-Gupta-Leimkuhler (GGL) formulation.
 *
 * For N coordinates and M constraints the integrated state is
 * \(z = [q_N, v_N, \lambda_M, \mu_M]^T\), so the @p nstates argument supplied
 * to the DAE integrators must equal \(2(N + M)\), and the initial conditions
 * must use this ordering.  The model is
 * \[
 * \dot q = v + G^T\mu, \qquad M(q)\dot v = F(x,q,v) - G^T\lambda,
 * \]
 * \[
 * 0 = G\dot q + \alpha g, \qquad
 * 0 = DG(q)[\dot q]v + GM^{-1}(F - G^T\lambda) + \beta Gv,
 * \]
 * where \(g(q)\) are the position constraints and
 * \(G = \partial g / \partial q\).  Every callback must accept perturbed,
 * off-manifold states, and the constraint callbacks must not depend on the
 * independent variable.
 */
typedef struct {
    /**
     * Required callback computing the applied force \(F(x,q,v)\), excluding
     * multiplier forces.  It is called with n = 2N and receives the
     * 2N-element vector [q, v], but must write only the N-element force
     * vector.
     */
    c_ode force;
    /**
     * Optional callback receiving the N-element coordinate vector q and
     * writing the N-by-N symmetric positive-definite physical mass matrix.
     * Set to NULL to use the identity matrix.  Set the mass_matrix_constant
     * setting to false whenever this matrix varies with x or q.
     */
    c_ode_mass_matrix mass;
    /**
     * Required callback receiving the N-element coordinate vector q and
     * writing the M position constraint residuals \(g(q)\).
     */
    c_dae_algebraic_constraints constraints;
    /**
     * Required callback receiving the N-element coordinate vector q and
     * writing the M-by-N constraint Jacobian \(G = \partial g / \partial q\)
     * in column-major order.  The rows of \(G\) must be linearly independent.
     */
    c_ode_jacobian gradient;
    /**
     * Required callback computing the M curvature terms \(DG(q)[\dot q]v\).
     * Return zeros for linear constraints.
     */
    c_dae_constraint_curvature curvature;
    /**
     * Nonnegative, finite position-constraint decay rate, in inverse
     * independent-variable units.  Zero disables position damping.
     */
    double alpha;
    /**
     * Nonnegative, finite velocity-constraint decay rate, in inverse
     * independent-variable units.  Zero disables velocity damping.
     */
    double beta;
    /** The number M of independent constraints, with 1 <= M <= N. */
    int constraint_count;
} c_ggl_container;

#ifdef __cplusplus
extern "C"{
#endif

/**
 * @brief Cleans up a c_ode_results type.
 * 
 * @param ptr The c_ode_results object to free.  It's the caller's
 * responsibility to nullify ptr.  This routine merely frees the
 * associated memory.
 */
void c_free_ode_results(c_ode_results ptr);

/**
 * @brief Gets the integrator performance statistics from a c_ode_results object.
 * 
 * @param ptr The c_ode_results object.
 * @param stats The c_ode_results_statistics object where the performance 
 * statistics will be written.
 * @return Returns C_NO_ERROR if successful; else C_NULL_POINTER_ERROR if ptr
 * is invalid or null.
 */
int c_get_ode_results_statistics(const c_ode_results ptr, c_ode_results_statistics *stats);

/**
 * @brief Gets the size of the matrix required to hold the integration results.
 * 
 * @param ptr The c_ode_results object.
 * @param rows The number of rows required to fit the results.
 * @param cols The number of columns required to fit the results.
 * @return Returns C_NO_ERROR if successful; else C_NULL_POINTER_ERROR if ptr
 * is invalid or null.
 */
int c_get_ode_results_size(const c_ode_results ptr, int *rows, int *cols);

/**
 * Returns the contents of the c_ode_results object as a matrix.
 * 
 * @param ptr The c_ode_results object.
 * @param m The number of rows in x.
 * @param n The number of columns in x.
 * @param x The M-by-N matrix where the results will be written.  The matrix is 
 * stored in column-major format with the first column containing the
 * independent variable values at which the integrator returned a solution 
 * point.  The remaining columns each contain a state variable.
 * @param ldx The matrix stride.
 * @return Returns C_NO_ERROR if successful; else, if ptr is invalid or null
 * a C_NULL_POINTER_ERROR is returned, if m is incorrectly sized (see 
 * c_get_ode_results_size) a C_MATRIX_ROW_SIZE_ERROR is returned, if n is
 * incorrectly sized (see c_get_ode_results_size) a C_MATRIX_COLUMN_SIZE_ERROR
 * is returned, and if ldx is not at least m a C_MATRIX_STRIDE_ERROR is returned.
 */
int c_get_ode_results(const c_ode_results ptr, int m, int n, double *x, int ldx);

/**
 * @brief Writes the contents of the c_ode_results object to a delimited text
 * file.
 * 
 * @param fname The filename.
 * @param ptr The c_ode_results object.
 * @param delimiter The delimiter character.
 * @return Returns C_NO_ERROR if successful; else C_NULL_POINTER_ERROR if ptr
 * is invalid or null.
 */
int c_write_ode_results(const char *fname, const c_ode_results ptr, char delimiter);

/**
 * @brief Fills the supplied c_ode_settings container with default settings.
 * 
 * @param settings The c_ode_settings container.
 */
void c_default_ode_settings(c_ode_settings *settings);

/**
 * @brief Utilizes the 3rd order Bogacki-Shampine Runge-Kutta integrator to 
 * solve a system of first order differential equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_runge_kutta_23(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_results *output);

/**
 * @brief Utilizes the 5th order Dormand-Prince Runge-Kutta integrator to 
 * solve a system of first order differential equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_runge_kutta_45(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_results *output);

/**
 * @brief Utilizes the 8th order Hairer Runge-Kutta integrator to 
 * solve a system of first order differential equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_runge_kutta_853(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_results *output);

/**
 * @brief Utilizes the 5th order Tsitouras Runge-Kutta integrator to 
 * solve a system of first order differential equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_tsitouras(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_results *output);

/**
 * @brief Utilizes a variable order Adams predictor-corrector integrator to 
 * solve a system of first order differential equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param mass A pointer to the routine for computing the mass matrix.  Pass
 * NULL if no routine is necessary.  This routine cannot handle singular mass
 * matrices.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_adams(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_mass_matrix mass, c_ode_results *output);

/**
 * @brief Utilizes a 4th order Rosenbrock integrator to solve a system of
 * first order differential equations or index-1 differential algebraic 
 * equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param jac A pointer to the routine for computing the Jacobian matrix.
 * Pass NULL if no routine is necessary and allowing the integrator to use 
 * a finite-difference estimate is acceptable.
 * @param mass A pointer to the routine for computing the mass matrix.  Pass
 * NULL if no routine is necessary.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_rosenbrock(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_jacobian jac, c_ode_mass_matrix mass, c_ode_results *output);

/**
 * @brief Utilizes a 4th order Kennedy-Carpenter integrator to solve a system of
 * first order differential equations or index-1 differential algebraic 
 * equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param jac A pointer to the routine for computing the Jacobian matrix.
 * Pass NULL if no routine is necessary and allowing the integrator to use 
 * a finite-difference estimate is acceptable.
 * @param mass A pointer to the routine for computing the mass matrix.  Pass
 * NULL if no routine is necessary.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_kennedy_carpenter_4(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_jacobian jac, c_ode_mass_matrix mass, c_ode_results *output);

/**
 * @brief Utilizes a 5th order Kennedy-Carpenter integrator to solve a system of
 * first order differential equations or index-1 differential algebraic 
 * equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param jac A pointer to the routine for computing the Jacobian matrix.
 * Pass NULL if no routine is necessary and allowing the integrator to use 
 * a finite-difference estimate is acceptable.
 * @param mass A pointer to the routine for computing the mass matrix.  Pass
 * NULL if no routine is necessary.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_kennedy_carpenter_5(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_jacobian jac, c_ode_mass_matrix mass, c_ode_results *output);

/**
 * @brief Utilizes a variable order BDF integrator to solve a system of
 * first order differential equations or index-1 differential algebraic 
 * equations.
 * 
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param fcn A pointer to the function containing the equations to integrate.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param jac A pointer to the routine for computing the Jacobian matrix.
 * Pass NULL if no routine is necessary and allowing the integrator to use 
 * a finite-difference estimate is acceptable.
 * @param mass A pointer to the routine for computing the mass matrix.  Pass
 * NULL if no routine is necessary.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_bdf(int nx, int nstates, c_ode fcn, const double *x, 
    const double *initcond, const c_ode_settings *settings, 
    c_ode_jacobian jac, c_ode_mass_matrix mass, c_ode_results *output);

/**
 * @brief Utilizes a 4th order Rosenbrock integrator to solve a system of
 * differential algebraic equations.
 * 
 * @param model The c_ggl_container containing the model and constraints.
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_dae_rosenbrock(const c_ggl_container *model, int nx, int nstates,
    const double *x, const double *initcond, const c_ode_settings *settings,
    c_ode_results *output);

/**
 * @brief Utilizes a 4th order Kennedy-Carpenter integrator to solve a system of
 * differential algebraic equations.
 * 
 * @param model The c_ggl_container containing the model and constraints.
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_dae_kennedy_carpenter_4(const c_ggl_container *model, int nx, int nstates,
    const double *x, const double *initcond, const c_ode_settings *settings,
    c_ode_results *output);

/**
 * @brief Utilizes a 5th order Kennedy-Carpenter integrator to solve a system of
 * differential algebraic equations.
 * 
 * @param model The c_ggl_container containing the model and constraints.
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_dae_kennedy_carpenter_5(const c_ggl_container *model, int nx, int nstates,
    const double *x, const double *initcond, const c_ode_settings *settings,
    c_ode_results *output);

/**
 * @brief Utilizes a variable order BDF integrator to solve a system of
 * differential algebraic equations.
 * 
 * @param model The c_ggl_container containing the model and constraints.
 * @param nx The number of solution points.  This value must be at least 2;
 * else, a C_INVALID_INPUT_ERROR is returned.
 * @param nstates The number of state equations.  This value must be at least 
 * 1; else, a C_INVALID_INPUT_ERROR is returned.
 * @param x An NX-element array containing the desired solution points.  If
 * only 2 points are defined, the integrator will return all successful 
 * intermediate integration points as well as the requested start and end 
 * points.  If more points are specified, the integrator will utilize dense
 * output and return the solution at the requested points via interpolation 
 * as required.
 * @param initcond An NSTATES-element array containing the initial conditions.
 * @param settings The integrator settings.  Use c_default_ode_settings to
 * utilize default settings.
 * @param output A pointer to the output c_ode_results handle.  Memory must be
 * freed by c_free_ode_results.
 * @return An error flag denoting the status of the operation.
 */
int c_dae_bdf(const c_ggl_container *model, int nx, int nstates,
    const double *x, const double *initcond, const c_ode_settings *settings,
    c_ode_results *output);

#ifdef __cplusplus
}
#endif
#endif
