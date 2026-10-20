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

#ifdef __cplusplus
}
#endif
#endif
