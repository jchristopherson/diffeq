!> @brief This module contains several ODE solvers and associated types.
module diffeq
    use iso_fortran_env
    use ferror
    implicit none
    private
    public :: ode
    public :: ode_jacobian
    public :: ode_mass_matrix
    public :: ode_solver
    public :: ode_integer_inquiry
    public :: ode_fixed_step
    public :: ode_container
    public :: ode_integrator
    public :: fixed_step_integrator
    public :: rk_fixed_integrator
    public :: rk4_fixed_integrator
    public :: fixed_multistep_integrator
    public :: adams_fixed_integerator
    public :: variable_step_integrator
    public :: variable_singlestep_integrator
    public :: rk_variable_integrator
    public :: dprk45_integrator
    public :: bsrk32_integrator
    public :: implicit_rk_variable_integrator
    public :: dirk_integrator
    public :: sdirk4_integrator
    public :: DIFFEQ_MEMORY_ALLOCATION_ERROR
    public :: DIFFEQ_NULL_POINTER_ERROR
    public :: DIFFEQ_MATRIX_SIZE_ERROR
    public :: DIFFEQ_ARRAY_SIZE_ERROR
    public :: DIFFEQ_INVALID_INPUT_ERROR
    public :: DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR
    public :: DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR
    public :: DIFFEQ_INVALID_OPERATION_ERROR

! ------------------------------------------------------------------------------
    integer(int32), parameter :: DIFFEQ_MEMORY_ALLOCATION_ERROR = 10000
    integer(int32), parameter :: DIFFEQ_NULL_POINTER_ERROR = 10001
    integer(int32), parameter :: DIFFEQ_MATRIX_SIZE_ERROR = 10002
    integer(int32), parameter :: DIFFEQ_ARRAY_SIZE_ERROR = 10003
    integer(int32), parameter :: DIFFEQ_INVALID_INPUT_ERROR = 10004
    integer(int32), parameter :: DIFFEQ_MISSING_ARGUMENT_ERROR = 10005
    integer(int32), parameter :: DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR = 10006
    integer(int32), parameter :: DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR = 10007
    integer(int32), parameter :: DIFFEQ_INVALID_OPERATION_ERROR = 10008

! ------------------------------------------------------------------------------
    !> @brief A container for the routine containing the ODEs to integrate.
    type ode_container
        ! A value determining if the mass matrix is state dependent such that it
        ! must be recomputed at each step.
        logical, private :: m_massDependent = .true.
        ! Jacobian calculation workspace array.
        real(real64), private, allocatable, dimension(:) :: m_jwork
        ! Finite difference step size.
        real(real64), private :: m_fdStep = sqrt(epsilon(1.0d0))
        !> @brief A pointer to the routine containing the ODEs to integrate.
        procedure(ode), pointer, public, nopass :: fcn => null()
        !> @brief A pointer to the routine containing the analytical Jacobian.
        !! If supplied, this routine is utilized; however, if null, a finite
        !! difference approximation is utilized.
        procedure(ode_jacobian), pointer, public, nopass :: &
            jacobian => null()
        !> @brief A pointer to the routine containing the mass matrix for the
        !! system.  If set to null (the default), an identity mass matrix will
        !! be assumed.
        procedure(ode_mass_matrix), pointer, public, nopass :: &
            mass_matrix => null()
    contains
        ! Use to allocate internal workspaces.  This routine only takes action
        ! if the workspace array(s) are not sized properly for the application.
        procedure, private :: allocate_workspace => oc_alloc_workspace
        !> @brief Gets the size of the step to use for the finite difference
        !! calculations used to estimate the Jacobian.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_finite_difference_step( &
        !!  class(ode_container) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref ode_container object.
        !! @return The step size.
        procedure, public :: get_finite_difference_step => oc_get_fd_step
        !> @brief Sets the size of the step to use for the finite difference
        !! calculations used to estimate the Jacobian.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_finite_difference_step( &
        !!  class(ode_container) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref ode_container object.
        !! @param[in] x The step size.
        procedure, public :: set_finite_difference_step => oc_set_fd_step
        !> @brief Gets a value determining if the mass matrix is state-dependent
        !! such that it requires updating at every integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function get_is_mass_matrix_dependent( &
        !!  class(ode_container) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref ode_container object.
        !! @return True if the mass matrix is state-dependent such that it 
        !!  requires updating at each integration step; else, false if the
        !!  mass matrix is not state-dependent and can be treated as constant
        !!  for all integration steps.
        procedure, public :: get_is_mass_matrix_dependent => &
            oc_get_is_mass_dependent
        !> @brief Sets a value determining if the mass matrix is state-dependent
        !! such that it requires updating at every integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_is_mass_matrix_dependent( &
        !!  class(ode_container) this, &
        !!  logical x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref ode_container object.
        !! @param[in] x True if the mass matrix is state-dependent such that it 
        !!  requires updating at each integration step; else, false if the
        !!  mass matrix is not state-dependent and can be treated as constant
        !!  for all integration steps.
        procedure, public :: set_is_mass_matrix_dependent => &
            oc_set_is_mass_dependent
        !> @brief Computes the Jacobian matrix for the system of ODEs.  If
        !! a routine is provided with an analytical Jacobian, the supplied
        !! routine is utilized; else, the Jacobian is estimated via a forward
        !! difference approximation.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine compute_jacobian( &
        !!  class(ode_container) this, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) jac(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref ode_container object.
        !! @param[in] x The current independent variable value.
        !! @param[in] y An N-element array containing the current dependent
        !!  variable values.
        !! @param[out] jac An N-by-N matrix where the Jacobian will be written.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling. Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined,
        !!      and the calculation is being performed by finite differences.
        !!  - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if @p jac is not N-by-N.
        procedure, public :: compute_jacobian => oc_jacobian
        !> @brief Evaluates the ODEs by evaluating the routine defined by
        !! @ref fcn.  This routine may also be overidden to provide custom
        !! functionallity.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine ode( &
        !!  class(ode_container) this, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) dydx(:) &
        !! )
        !!
        !! @param[in] this The @ref ode_container object.
        !! @param[in] x The current value of the independent variables.
        !! @param[in] y An N-element array containing the current values of the
        !!  dependent variables.
        !! @param[out] dydx An N-element array where the output of each of the
        !!  N ODEs will be written.
        procedure, public :: ode => oc_ode_fcn
    end type

    ! diffeq_ode_container.f90
    interface
        subroutine ode(x, y, dydx)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        subroutine ode_jacobian(x, y, jac)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        subroutine ode_mass_matrix(x, y, m)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: m
        end subroutine

        pure module function oc_get_fd_step(this) result(rst)
            class(ode_container), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine oc_set_fd_step(this, x)
            class(ode_container), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function oc_get_is_mass_dependent(this) result(rst)
            class(ode_container), intent(in) :: this
            logical :: rst
        end function

        module subroutine oc_set_is_mass_dependent(this, x)
            class(ode_container), intent(inout) :: this
            logical :: x
        end subroutine

        module subroutine oc_alloc_workspace(this, ndof, err)
            class(ode_container), intent(inout) :: this
            integer(int32), intent(in) :: ndof
            class(errors), intent(inout) :: err
        end subroutine

        module subroutine oc_jacobian(this, x, y, jac, err)
            class(ode_container), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine oc_ode_fcn(this, x, y, dydx)
            class(ode_container), intent(in) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief The most basic ODE integrator object capable of integrating
    !! systems of ODE's.
    type, abstract :: ode_integrator
    contains
        !> @brief Solves the supplied system of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! allocatable real(real64)(:,:) function solve( &
        !!  class(ode_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) x(:), &
        !!  real(real64) iv(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref ode_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] x An array, of at least 2 values, defining at a minimum
        !!  the starting and ending values of the independent variable 
        !!  integration range.  If more than two values are specified, the
        !!  integration results will be returned at the supplied values.
        !! @param[in] An array containing the initial values for each ODE.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined.
        !!  - DIFFEQ_INVALID_INPUT_ERROR: Occurs if there are less than 2 values
        !!      given in the independent variable array @p x.
        !!
        !! @return An M-by-N matrix where M is the number of solution points, 
        !!  and N is the number of ODEs plus 1.  The first column contains
        !!  the values of the independent variable at which the results were
        !!  computed.  The remaining columns contain the integration results
        !!  for each ODE.
        procedure(ode_solver), public, pass, deferred :: solve
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref ode_integrator object.
        !! @return The order of the integrator.
        procedure(ode_integer_inquiry), public, pass, deferred :: get_order
    end type

    interface
        function ode_solver(this, sys, x, iv, err) result(rst)
            use iso_fortran_env
            use ferror
            import ode_integrator
            import ode_container
            class(ode_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function

        pure function ode_integer_inquiry(this) result(rst)
            use iso_fortran_env
            import ode_integrator
            class(ode_integrator), intent(in) :: this
            integer(int32) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a fixed-step integrator.
    type, abstract, extends(ode_integrator) :: fixed_step_integrator
    contains
        !> @brief Solves the supplied system of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! allocatable real(real64)(:,:) function solve( &
        !!  class(fixed_step_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) x(:), &
        !!  real(real64) iv(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref fixed_step_integrator object.
        !! @param[in,out] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] x An array, of at least 2 values, defining at a minimum
        !!  the starting and ending values of the independent variable 
        !!  integration range.  If more than two values are specified, the
        !!  integration results will be returned at the supplied values.
        !! @param[in] An array containing the initial values for each ODE.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        !!
        !! @return An M-by-N matrix where M is the number of solution points, 
        !!  and N is the number of ODEs plus 1.  The first column contains
        !!  the values of the independent variable at which the results were
        !!  computed.  The remaining columns contain the integration results
        !!  for each ODE.
        procedure, public :: solve => fsi_solver
        !> @brief Takes one integration step of a predetermined size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine step( &
        !!  class(fixed_step_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional real(real64) fprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref fixed_step_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[out] fprev An M-by-NEQN array where the previous M function
        !!  values are written.  This is typically useful for multi-step 
        !!  methods.  In single-step methods this parameter is not used.
        !! @param[in,out] An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure(ode_fixed_step), public, pass, deferred :: step
    end type

    ! diffeq_fs_integrator.f90
    interface
        subroutine ode_fixed_step(this, sys, h, x, y, yn, xprev, yprev, fprev, &
            err)
            use iso_fortran_env 
            use ferror   
            import fixed_step_integrator
            import ode_container
            class(fixed_step_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module function fsi_solver(this, sys, x, iv, err) result(rst)
            class(fixed_step_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines an explicit, Runge-Kutta fixed-step integrator.
    type, abstract, extends(fixed_step_integrator) :: rk_fixed_integrator
        ! Workspace matrix
        real(real64), private, allocatable, dimension(:,:) :: m_work
    contains
        ! Use to allocate internal workspaces.  This routine only takes action
        ! if the workspace array(s) are not sized properly for the application.
        procedure, private :: allocate_workspace => rkf_alloc_workspace
        !> @brief Takes one integration step of a predetermined size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine step( &
        !!  class(rk_fixed_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_fixed_integrator object.
        !! @param[in,out] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[in,out] An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: step => rkf_step
        !> @brief Gets the requested method factor from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_method_factor( &
        !!  class(rk_fixed_integrator) this, &
        !!  integer(int32), intent(in) i, &
        !!  integer(int32), intent(in) j &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk_fixed_integrator object.
        !! @param[in] i The row index of the parameter from the Butcher tableau.
        !! @param[in] j The column index of the parameter from the Butcher 
        !!  tableau.
        !! @return The requested parameter.
        procedure(rkf_get_matrix_parameter), deferred, public :: &
            get_method_factor
        !> @brief Gets the requested quadrature weight from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_quadrature_weight( &
        !!  class(rk_fixed_integrator) this, &
        !!  integer(int32), intent(in) i &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk_fixed_integrator object.
        !! @param[in] i The index of the parameter from the Butcher tableau.
        !! @return The requested parameter.
        procedure(rkf_get_array_parameter), deferred, public :: &
            get_quadrature_weight
        !> @brief Gets the requested position factor from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_position_factor( &
        !!  class(rk_fixed_integrator) this, &
        !!  integer(int32), intent(in) i &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk_fixed_integrator object.
        !! @param[in] i The index of the parameter from the Butcher tableau.
        !! @return The requested parameter.
        procedure(rkf_get_array_parameter), deferred, public :: &
            get_position_factor
    end type

    interface
        pure function rkf_get_matrix_parameter(this, i, j) result(rst)
            use iso_fortran_env
            import rk_fixed_integrator
            class(rk_fixed_integrator), intent(in) :: this
            integer(int32), intent(in) :: i, j
            real(real64) :: rst
        end function

        pure function rkf_get_array_parameter(this, i) result(rst)
            use iso_fortran_env
            import rk_fixed_integrator
            class(rk_fixed_integrator), intent(in) :: this
            integer(int32), intent(in) :: i
            real(real64) :: rst
        end function

        module subroutine rkf_alloc_workspace(this, n, neqn, err)
            class(rk_fixed_integrator), intent(inout) :: this
            integer(int32), intent(in) :: n, neqn
            class(errors), intent(inout) :: err
        end subroutine

        module subroutine rkf_step(this, sys, h, x, y, yn, xprev, yprev, &
            fprev, err)
            class(rk_fixed_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the explicit, 4th order, Runge-Kutta fixed-step 
    !! integrator.
    !!
    !! @par Example
    !! The following example solves the first-order differential equation
    !! \f$ \frac{dy}{dx} + y \sin^{2}x = 0 \f$ where the solution is 
    !! \f$ y = C \exp \left( \frac{1}{4} \sin(2 x) - \frac{x}{2} \right) \f$.
    !!
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq
    !!     use diffeq_models
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: npts = 10000
    !!     real(real64), parameter :: h = 1.0d-4
    !!
    !!     ! Local Variables
    !!     type(rk4_fixed_integrator) :: integrator
    !!     type(ode_container) :: mdl
    !!     integer(int32) :: i
    !!     real(real64) :: x(npts), analytical(npts), numerical(npts, 2)
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd1, pd2
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Define the values of x at which the solution is to be computed
    !!     x = (/ (i * h, i = 0, npts - 1) /)
    !!
    !!     ! Define the model
    !!     mdl%fcn => test_1dof_1
    !!
    !!     ! Compute the solution, both numerical and analytical
    !!     numerical = integrator%solve(mdl, x, [2.0d0])
    !!     analytical = test_1dof_solution_1(x)
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     lgnd => plt%get_legend()
    !!     call xAxis%set_title("x")
    !!     call yAxis%set_title("y(x)")
    !!     call lgnd%set_is_visible(.true.)
    !!     call lgnd%set_horizontal_position(LEGEND_LEFT)
    !!     call lgnd%set_vertical_position(LEGEND_BOTTOM)
    !!     call plt%set_title("y' + y sin^2 x = 0 ")
    !!
    !!     call pd1%define_data(x, analytical)
    !!     call pd1%set_name("Analytical")
    !!     call pd1%set_line_width(2.0)
    !!     call plt%push(pd1)
    !!
    !!     call pd2%define_data(numerical(:,1), numerical(:,2))
    !!     call pd2%set_name("Numerical")
    !!     call pd2%set_line_width(4.0)
    !!     call pd2%set_line_style(LINE_DASHED)
    !!     call plt%push(pd2)
    !!
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The ODE routine was stored in a seperate module; however, here is the
    !! code for the ODE routine.
    !! @code{.f90}
    !! subroutine test_1dof_1(x, y, dydx)
    !!     real(real64), intent(in) :: x, y(:)
    !!     real(real64), intent(out) :: dydx(:)
    !!     dydx(1) = -y(1) * sin(x)**2
    !! end subroutine
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html rk4_example_1.png
    type, extends(rk_fixed_integrator) :: rk4_fixed_integrator
    contains
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order(class(rk4_fixed_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref rk4_fixed_integrator object.
        !! @return The order of the integrator.
        procedure, public :: get_order => rk4_get_order
        !> @brief Gets the requested method factor from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_method_factor( &
        !!  class(rk4_fixed_integrator) this, &
        !!  integer(int32), intent(in) i, &
        !!  integer(int32), intent(in) j &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk4_fixed_integrator object.
        !! @param[in] i The row index of the parameter from the Butcher tableau.
        !! @param[in] j The column index of the parameter from the Butcher 
        !!  tableau.
        !! @return The requested parameter.
        procedure, public :: get_method_factor => rk4_get_a_param
        !> @brief Gets the requested quadrature weight from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_quadrature_weight( &
        !!  class(rk4_fixed_integrator) this, &
        !!  integer(int32), intent(in) i &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk4_fixed_integrator object.
        !! @param[in] i The index of the parameter from the Butcher tableau.
        !! @return The requested parameter.
        procedure, public :: get_quadrature_weight => rk4_get_b_param
        !> @brief Gets the requested position factor from the Butcher tableau.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_position_factor( &
        !!  class(rk4_fixed_integrator) this, &
        !!  integer(int32), intent(in) i &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk4_fixed_integrator object.
        !! @param[in] i The index of the parameter from the Butcher tableau.
        !! @return The requested parameter.
        procedure, public :: get_position_factor => rk4_get_c_param
    end type

    ! diffeq_rk4.f90
    interface
        pure module function rk4_get_order(this) result(rst)
            class(rk4_fixed_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function rk4_get_a_param(this, i, j) result(rst)
            class(rk4_fixed_integrator), intent(in) :: this
            integer(int32), intent(in) :: i, j
            real(real64) :: rst
        end function

        pure module function rk4_get_b_param(this, i) result(rst)
            class(rk4_fixed_integrator), intent(in) :: this
            integer(int32), intent(in) :: i
            real(real64) :: rst
        end function

        pure module function rk4_get_c_param(this, i) result(rst)
            class(rk4_fixed_integrator), intent(in) :: this
            integer(int32), intent(in) :: i
            real(real64) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a fixed step-size, multi-step integrator.
    type, abstract, extends(fixed_step_integrator) :: fixed_multistep_integrator
        ! An NEQN-by-ORDER storage matrix for ODE outputs.
        real(real64), allocatable, dimension(:,:) :: m_buffer
    contains
        ! Use to allocate internal workspaces.  This routine only takes action
        ! if the workspace array(s) are not sized properly for the application.
        procedure, private :: allocate_workspace => fms_alloc_workspace
        !> @brief Solves the supplied system of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! allocatable real(real64)(:,:) function solve( &
        !!  class(fixed_multistep_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) x(:), &
        !!  real(real64) iv(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref fixed_multistep_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] x An array, of at least 2 values, defining at a minimum
        !!  the starting and ending values of the independent variable 
        !!  integration range.  If more than two values are specified, the
        !!  integration results will be returned at the supplied values.
        !! @param[in] An array containing the initial values for each ODE.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE function is defined.
        !!  - DIFFEQ_INVALID_INPUT_ERROR: Occurs if there are less than 2 values
        !!      given in the independent variable array @p x.
        !!
        !! @return An M-by-N matrix where M is the number of solution points, 
        !!  and N is the number of ODEs plus 1.  The first column contains
        !!  the values of the independent variable at which the results were
        !!  computed.  The remaining columns contain the integration results
        !!  for each ODE.
        procedure, public :: solve => fms_solver
    end type

    ! diffeq_multistep_fixed.f90
    interface
        module subroutine fms_alloc_workspace(this, neqn, err)
            class(fixed_multistep_integrator), intent(inout) :: this
            integer(int32), intent(in) :: neqn
            class(errors), intent(inout) :: err
        end subroutine
        
        module function fms_solver(this, sys, x, iv, err) result(rst)
            class(fixed_multistep_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a fixed-step, 4th-order, Adams-Bashforth-Moulton PECE 
    !! integrator.
    !!
    !! @par Example
    !! The following example solves the Van der Pol equation \f$ 
    !! \frac{d^2y}{dx^2} - \mu \left( 1 - y^2 \right) \frac{dy}{dx} + y = 0 \f$
    !! using both this Adams method integrator, and the 4th order Runge-Kutta
    !! integrator @ref rk4_fixed_integrator.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq
    !!     use diffeq_models
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: npts = 25000
    !!     real(real64), parameter :: h = 1.0d-3
    !!
    !!     ! Local Variables
    !!     type(adams_fixed_integerator) :: integrator
    !!     type(rk4_fixed_integrator) :: comparison
    !!     type(ode_container) :: mdl
    !!     integer(int32) :: i
    !!     real(real64) :: x(npts), sol(npts, 3), csol(npts, 3)
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd1, pd2
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Define the values of x at which the solution is to be computed
    !!     x = (/ (i * h, i = 0, npts - 1) /)
    !!
    !!     ! Define the model
    !!     mdl%fcn => vanderpol
    !!
    !!     ! Compute the solution with both integrators
    !!     sol = integrator%solve(mdl, x, [2.0d0, 0.0d0])
    !!     csol = comparison%solve(mdl, x, [2.0d0, 0.0d0])
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     lgnd => plt%get_legend()
    !!     call xAxis%set_title("x")
    !!     call yAxis%set_title("y(x)")
    !!     call lgnd%set_is_visible(.true.)
    !!     call lgnd%set_horizontal_position(LEGEND_CENTER)
    !!     call lgnd%set_vertical_position(LEGEND_BOTTOM)
    !!
    !!     call pd1%define_data(sol(:,1), sol(:,2))
    !!     call pd1%set_name("Adams")
    !!     call pd1%set_line_width(2.0)
    !!     call plt%push(pd1)
    !!
    !!     call pd2%define_data(csol(:,1), csol(:,2))
    !!     call pd2%set_name("Runge Kutta")
    !!     call pd2%set_line_width(4.0)
    !!     call pd2%set_line_style(LINE_DASHED)
    !!     call plt%push(pd2)
    !!
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The ODE routine was stored in a seperate module; however, here is the
    !! code for the ODE routine.
    !! @code{.f90}
    !! subroutine vanderpol(x, y, dydx)
    !!     ! Arguments
    !!     real(real64), intent(in) :: x, y(:)
    !!     real(real64), intent(out) :: dydx(:)
    !!
    !!     ! Model Constants
    !!     real(real64), parameter :: mu = 5.0d0
    !!
    !!     ! Equations
    !!     dydx(1) = y(2)
    !!     dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
    !! end subroutine
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html vanderpol_fixed_adams_vs_runge_kutta_example_1.png
    type, extends(fixed_multistep_integrator) :: adams_fixed_integerator
    contains
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order(class(adams_fixed_integerator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref adams_fixed_integerator object.
        !! @return The order of the integrator.
        procedure, public :: get_order => afi_get_order
        !> @brief Takes one integration step of a predetermined size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine step( &
        !!  class(adams_fixed_integerator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref adams_fixed_integerator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !! @param[out] fprev An M-by-NEQN array where the previous M function
        !!  values are written.  This is typically useful for multi-step 
        !!  methods.  In single-step methods this parameter is not used.
        !! @param[in,out] An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: step => afi_step
    end type

    ! diffeq_abmf.f90
    interface
        pure module function afi_get_order(this) result(rst)
            class(adams_fixed_integerator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine afi_step(this, sys, h, x, y, yn, xprev, yprev, &
            fprev, err)
            class(adams_fixed_integerator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface





! ******************************************************************************
! ******************************************************************************
! VARIABLE STEP INTEGRATORS
! ******************************************************************************
! ******************************************************************************





    !> @brief Defines a variable-step integrator.
    type, abstract, extends(ode_integrator) :: variable_step_integrator
        real(real64), private :: m_safetyfactor = 0.9d0
        real(real64), private :: m_maxstep = huge(1.0d0)
        real(real64), private :: m_minstep = 1.0d2 * epsilon(1.0d0)
        integer(int32), private :: m_maxitercount = 1000
        integer(int32), private :: m_maxstepcount = 1000000
        real(real64), private :: m_stepSize = 1.0d0
        real(real64), private :: m_nextStep = 1.0d0
        real(real64), private :: m_enormPrev = 1.0d0
        logical, private :: m_respectXMax = .true.
        ! Solution buffer
        real(real64), private, allocatable, dimension(:,:) :: m_buffer
        integer(int32) :: m_bufferCount = 0
        ! Workspaces
        real(real64), allocatable, dimension(:) :: m_ework      ! NEQN element
        ! Tolerances
        real(real64), allocatable, dimension(:) :: m_rtol       ! NEQN element
        real(real64), allocatable, dimension(:) :: m_atol       ! NEQN element
    contains
        !> @brief Initializes the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine initialize( &
        !!  class(variable_step_integrator) this, &
        !!  integer(int32) neqn, &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] neqn The number of equations being integrated.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        procedure, public :: initialize => vsi_alloc_workspace
        !> @brief Attempts a single integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine attempt_step( &
        !!  class(variable_step_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) en(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional real(real64) fprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[out] en An N-element array where the values of the error
        !!  estimates will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] fprev An M-by-NEQN array containing the previous M arrays
        !!  of function (derivative) values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure(variable_step_attempt), deferred, public :: attempt_step
        !> @brief Perform necessary actions on completion of a successful step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine on_successful_step( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The new value of the independent variable.
        !! @param[in] y An N-element array containing the values of the 
        !!  dependent variables at @p x.
        !! @param[in] yn An N-element array containing the values of the
        !!  dependent variables at @p xn.
        procedure(variable_step_action), deferred, public :: on_successful_step
        !> @brief Gets a safety factor used to limit the predicted step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_safety_factor( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The parameter value.
        procedure, public :: get_safety_factor => vsi_get_safety_factor
        !> @brief Sets a safety factor used to limit the predicted step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_safety_factor( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_safety_factor => vsi_set_safety_factor
        !> @brief Gets the maximum allowed step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_max_step_size( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The parameter value.
        procedure, public :: get_max_step_size => vsi_get_max_step
        !> @brief Sets the maximum allowed step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_max_step_size => vsi_set_max_step
        !> @brief Gets the minimum allowed step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_min_step_size( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The parameter value.
        procedure, public :: get_min_step_size => vsi_get_min_step
        !> @brief Sets the minimum allowed step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_min_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_min_step_size => vsi_set_min_step
        !> @brief Gets the maximum number of iterations per step allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_max_per_step_iteration_count( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The parameter value.
        procedure, public :: get_max_per_step_iteration_count => &
            vsi_get_max_iter_count
        !> @brief Sets the maximum number of iterations per step allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_per_step_iteration_count( &
        !!  class(variable_step_integrator) this, &
        !!  integer(int32) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_max_per_step_iteration_count => &
            vsi_set_max_iter_count
        !> @brief Gets the maximum number of integration steps allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_max_integration_step_count( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The parameter value.
        procedure, public :: get_max_integration_step_count => & 
            vsi_get_max_step_count
        !> @brief Sets the maximum number of integration steps allowed.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_integration_step_count( &
        !!  class(variable_step_integrator) this, &
        !!  integer(int32) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x The parameter value.
        procedure, public :: set_max_integration_step_count => &
            vsi_set_max_step_count
        !> @brief Computes the next step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function compute_next_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) hn, &
        !!  real(real64) en, &
        !!  real(real64) enm1 &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] hn The current step size.
        !! @param[in] en The norm of the error for the current step size.
        !! @param[in] enm1 The norm of the error from the previous step size.
        !!
        !! @return The new step size.
        procedure(next_step_size_calculator), deferred, public :: &
            compute_next_step_size
        !> @brief Buffers a results set.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine buffer_results( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x The independent variable value.
        !! @param[in] y An N-element array containing the solution values.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if @p y is not compatible with
        !!      the buffer size.
        procedure, public :: buffer_results => vsi_append_to_buffer
        !> @brief Gets the number of entries into the solution buffer.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_buffer_size( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The number of buffer entries.
        procedure, public :: get_buffer_size => vsi_get_buffer_count
        !> @brief Clears the results buffer.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine clear_buffer(class(variable_step_integrator) this)
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        procedure, public :: clear_buffer => vsi_clear_buffer
        !> @brief Gets the current step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_step_size( &
        !!  class(variable_step_integrator) this &
        !! )
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The step size.
        procedure, public :: get_step_size => vsi_get_step_size
        !> @brief Sets the current step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x &
        !! )
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x The step size.
        procedure, public :: set_step_size => vsi_set_step_size
        !> @brief Gets the next step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_step_size( &
        !!  class(variable_step_integrator) this &
        !! )
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The step size.
        procedure, public :: get_next_step_size => vsi_get_next_step_size
        !> @brief Sets the next step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) x &
        !! )
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x The step size.
        procedure, public :: set_next_step_size => vsi_set_next_step_size
        !> @brief Gets a value determining if the integrator should respect a
        !! hard limit in the independent variable range.  If false, the 
        !! integrator may step pass the limit.
        !! 
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function get_respect_x_max( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return True if the integrator should respect the limiting value of
        !!  the independent variable; else, false.
        procedure, public :: get_respect_x_max => vsi_get_respect_xmax
        !> @brief Sets a value determining if the integrator should respect a
        !! hard limit in the independent variable range.  If false, the 
        !! integrator may step pass the limit.
        !! 
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_respect_x_max( &
        !!  class(variable_step_integrator) this, &
        !!  logical x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] x True if the integrator should respect the limiting 
        !!  value of the independent variable; else, false.
        procedure, public :: set_respect_x_max => vsi_set_respect_xmax
        !> @brief Takes one integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine step( &
        !!  class(variable_step_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) x, &
        !!  real(real64) xmax, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional real(real64) fprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_step_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xmax The upper integration limit.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[out] fprev An M-by-NEQN array where the previous M function
        !!  values are written.  This is typically useful for multi-step 
        !!  methods.  In single-step methods this parameter is not used.
        !! @param[in,out] An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if @p yn is not the same size
        !!      as @p y.
        !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size becomes
        !!      too small.
        !!  - DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR: Occurs if the iteration
        !!      count is exceeded for a single step.
        procedure, public :: step => vsi_step
        !> @brief Computes an estimate to the first step size based upon the
        !! initial function values.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function estimate_first_step_size( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) xo, &
        !!  real(real64) xf, &
        !!  real(real64) yo(:), &
        !!  real(real64) fo(:) &
        !! )
        !! @endcode
        !!
        !! @param[in] this The variable_step_integrator object.
        !! @param[in] xo The initial value of the independent variable.
        !! @param[in] xf The final value of the independent variable.
        !! @param[in] yo An N-element array containing the initial values.
        !! @param[in] fo An N-element array containing the initial function 
        !!  values.
        !!
        !! @return An estimate on the initial step size.
        procedure, public :: estimate_first_step_size => vsi_estimate_first_step
        !> @brief Provides interpolation between integration points allowing
        !! for dense output.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine interpolate( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) xprev, &
        !!  real(real64) xnew, &
        !!  real(real64) x, &
        !!  real(real64) y(:),
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The variable_step_integrator object.
        !! @param[in] xprev The previous value of the independent variable.
        !! @param[in] xnew The updated value of the independent variable.
        !! @param[in] x The value at which to perform the interpolation.
        !! @param[out] y An N-element array containing the interpolated 
        !!  values for each equation.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure(variable_step_interpolation), public, deferred :: interpolate
        !> @brief Gets the default relative error tolerance.
        !! 
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_default_relative_tolerance( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The tolerance value.
        procedure, public :: get_default_relative_tolerance => &
            vsi_get_default_rel_tol
        !> @brief Gets the default absolute error tolerance.
        !! 
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_default_absolute_tolerance( &
        !!  class(variable_step_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref variable_step_integrator object.
        !! @return The tolerance value.
        procedure, public :: get_default_absolute_tolerance => &
            vsi_get_default_abs_tol
    end type

    ! diffeq_vs_integrator.f90
    interface
        subroutine variable_step_attempt(this, sys, h, x, y, yn, en, xprev, &
            yprev, fprev, err)
            use iso_fortran_env
            use ferror
            import variable_step_integrator
            import ode_container
            class(variable_step_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn, en
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine

        subroutine variable_step_action(this, x, xn, y, yn)
            use iso_fortran_env
            import variable_step_integrator
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
        end subroutine

        subroutine variable_step_interpolation(this, xprev, xnew, x, y, err)
            use iso_fortran_env
            use ferror
            import variable_step_integrator
            class(variable_step_integrator), intent(in) :: this
            real(real64), intent(in) :: xprev, xnew, x
            real(real64), intent(out), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine

        function next_step_size_calculator(this, hn, en, enm1) result(rst)
            use iso_fortran_env
            import variable_step_integrator
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: hn, en, enm1
            real(real64) :: rst
        end function

        pure module function vsi_get_safety_factor(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_safety_factor(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function vsi_get_max_step(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_max_step(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function vsi_get_min_step(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_min_step(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function vsi_get_max_iter_count(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine vsi_set_max_iter_count(this, x)
            class(variable_step_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine

        pure module function vsi_get_max_step_count(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine vsi_set_max_step_count(this, x)
            class(variable_step_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine

        pure module function estimate_error_1(y, ys, atol, rtol) result(rst)
            real(real64), intent(in), dimension(:) :: y, ys, atol, rtol
            real(real64) :: rst
        end function

        pure module function estimate_error_2(y, ys, atol, rtol) result(rst)
            real(real64), intent(in), dimension(:) :: y, ys
            real(real64), intent(in) :: atol, rtol
            real(real64) :: rst
        end function

        module subroutine vsi_append_to_buffer(this, x, y, err)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in) :: y(:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function vsi_get_buffer_count(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine vsi_clear_buffer(this)
            class(variable_step_integrator), intent(inout) :: this
        end subroutine

        pure module function vsi_get_step_size(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_step_size(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function vsi_get_next_step_size(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_next_step_size(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        module subroutine vsi_alloc_workspace(this, neqn, err)
            class(variable_step_integrator), intent(inout) :: this
            integer(int32), intent(in) :: neqn
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function vsi_get_respect_xmax(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            logical :: rst
        end function

        module subroutine vsi_set_respect_xmax(this, x)
            class(variable_step_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        module subroutine vsi_step(this, sys, x, xmax, y, yn, xprev, yprev, &
            fprev, err)
            class(variable_step_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: x, xmax
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function vsi_estimate_first_step(this, xo, xf, yo, fo) &
        result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64), intent(in) :: xo, xf
            real(real64), intent(in), dimension(:) :: yo, fo
            real(real64) :: rst
        end function

        pure module function vsi_get_default_rel_tol(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        pure module function vsi_get_default_abs_tol(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a variable-step, single-stage integrator.
    type, abstract, extends(variable_step_integrator) :: &
        variable_singlestep_integrator
    contains
        !> @brief Solves the supplied system of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! allocatable real(real64)(:,:) function solve( &
        !!  class(variable_singlestep_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) x(:), &
        !!  real(real64) iv(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref variable_singlestep_integrator object.
        !! @param[in,out] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] x An array, of at least 2 values, defining at a minimum
        !!  the starting and ending values of the independent variable 
        !!  integration range.  If more than two values are specified, the
        !!  integration results will be returned at the supplied values.
        !! @param[in] An array containing the initial values for each ODE.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        !!  - DIFFEQ_ARRAY_SIZE_ERROR: Occurs if @p x has less than 2 elements.
        !!  - DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR: Occurs if the step size becomes
        !!      too small.
        !!  - DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR: Occurs if the iteration
        !!      count is exceeded for a single step.
        !!  - DIFFEQ_NULL_POINTER_ERROR: Occurs if no ODE routine is defined.
        !!  - DIFFEQ_INVALID_INPUT_ERROR: Occurs if max(@p x) - min(@p x) = 0.
        !!
        !! @return An M-by-N matrix where M is the number of solution points, 
        !!  and N is the number of ODEs plus 1.  The first column contains
        !!  the values of the independent variable at which the results were
        !!  computed.  The remaining columns contain the integration results
        !!  for each ODE.
        procedure, public :: solve => vssi_solve
        procedure, private :: solve_driver => vssi_solve_driver
        procedure, private :: dense_solve_driver => vssi_dense_solve_driver
    end type

    ! diffeq_vssi.f90
    interface
        module function vssi_solve(this, sys, x, iv, err) result(rst)
            class(variable_singlestep_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function

        module function vssi_solve_driver(this, sys, x, iv, err) result(rst)
            class(variable_singlestep_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout) :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function

        module function vssi_dense_solve_driver(this, sys, x, iv, err) &
            result(rst)
            class(variable_singlestep_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout) :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a variable-step, Runge-Kutta integrator.
    type, abstract, extends(variable_singlestep_integrator) :: &
        rk_variable_integrator
        ! Workspace matrix (NEQN -by- STAGE COUNT)
        real(real64), private, allocatable, dimension(:,:) :: m_work
        ! Workspace array (NEQN)
        real(real64), private, allocatable, dimension(:) :: m_ywork
        ! A flag determining if this is the first accepted step (use for FSAL)
        logical :: m_firstStep = .true.
        ! Step-size PI control parameters
        real(real64), private :: m_alpha = 0.7d0
        real(real64), private :: m_beta = 0.4d-1
        !> @brief The NSTAGE-by-NSTAGE method factor matrix A.
        real(real64), public, allocatable, dimension(:,:) :: a
        !> @brief The NSTAGE-element quadrature weight array B.
        real(real64), public, allocatable, dimension(:) :: b
        !> @brief The NSTAGE-element position factor array C.
        real(real64), public, allocatable, dimension(:) :: c
        !> @brief The NSTAGE-element error estimate factor array E.
        real(real64), public, allocatable, dimension(:) :: e
    contains
        !> @brief Initializes the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine initialize( &
        !!  class(rk_variable_integrator) this, &
        !!  integer(int32) neqn, &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] neqn The number of equations being integrated.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        procedure, public :: initialize => rkv_alloc_workspace
        ! !> @brief Gets the requested method factor from the Butcher tableau.
        ! !!
        ! !! @par Syntax
        ! !! @code{.f90}
        ! !! real(real64) pure function get_method_factor( &
        ! !!  class(rk_variable_integrator) this, &
        ! !!  integer(int32), intent(in) i, &
        ! !!  integer(int32), intent(in) j &
        ! !! )
        ! !! @endcode
        ! !!
        ! !! @param[in] this The @ref rk_variable_integrator object.
        ! !! @param[in] i The row index of the parameter from the Butcher tableau.
        ! !! @param[in] j The column index of the parameter from the Butcher 
        ! !!  tableau.
        ! !! @return The requested parameter.
        ! procedure(rkv_get_matrix_parameter), deferred, public :: &
        !     get_method_factor
        ! !> @brief Gets the requested quadrature weight from the Butcher tableau.
        ! !!
        ! !! @par Syntax
        ! !! @code{.f90}
        ! !! real(real64) pure function get_quadrature_weight( &
        ! !!  class(rk_variable_integrator) this, &
        ! !!  integer(int32), intent(in) i &
        ! !! )
        ! !! @endcode
        ! !!
        ! !! @param[in] this The @ref rk_variable_integrator object.
        ! !! @param[in] i The index of the parameter from the Butcher tableau.
        ! !! @return The requested parameter.
        ! procedure(rkv_get_array_parameter), deferred, public :: &
        !     get_quadrature_weight
        ! !> @brief Gets the requested error coefficient from the Butcher tableau.
        ! !!
        ! !! @par Syntax
        ! !! @code{.f90}
        ! !! real(real64) pure function get_error_factor( &
        ! !!  class(rk_variable_integrator) this, &
        ! !!  integer(int32), intent(in) i &
        ! !! )
        ! !! @endcode
        ! !!
        ! !! @param[in] this The @ref rk_variable_integrator object.
        ! !! @param[in] i The index of the parameter from the Butcher tableau.
        ! !! @return The requested parameter.
        ! procedure(rkv_get_array_parameter), deferred, public :: &
        !     get_error_factor
        ! !> @brief Gets the requested position factor from the Butcher tableau.
        ! !!
        ! !! @par Syntax
        ! !! @code{.f90}
        ! !! real(real64) pure function get_position_factor( &
        ! !!  class(rk_variable_integrator) this, &
        ! !!  integer(int32), intent(in) i &
        ! !! )
        ! !! @endcode
        ! !!
        ! !! @param[in] this The @ref rk_variable_integrator object.
        ! !! @param[in] i The index of the parameter from the Butcher tableau.
        ! !! @return The requested parameter.
        ! procedure(rkv_get_array_parameter), deferred, public :: &
        !     get_position_factor
        !> @brief Determines if the integrator is an FSAL (first same as last)
        !! integrator (e.g. the 4th/5th order Dormand-Prince integrator).
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function is_fsal(class(rk_variable_integrator) this)
        !! @endcode
        !! 
        !! @param[in] this The @ref rk_variable_integrator object.
        !! @return Returns true if the integrator is an FSAL type; else,
        !!  returns false.
        procedure(rkv_get_boolean_parameter), deferred, public :: is_fsal
        !> @brief Gets the number of stages used by the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_stage_count( &
        !!  class(rk_variable_integrator) this &
        !! )
        !! @endcode
        !! 
        !! @param[in] this The @ref rk_variable_integrator object.
        !! @return The number of stages.
        procedure(rkv_get_integer_parameter), deferred, public :: &
            get_stage_count
        !> @brief Defines (initializes) the model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_model(class(rkv_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref rkv_integrator object.
        procedure(rkv_action), deferred, public :: define_model
        !> @brief Resets the integrator to its initial state.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine reset(class(rk_variable_integrator) this)
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        procedure, public :: reset => rkv_reset
        !> @brief Attempts a single integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine attempt_step( &
        !!  class(rk_variable_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) en(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional real(real64) fprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[out] en An N-element array where the values of the error
        !!  estimates will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] fprev An M-by-NEQN array containing the previous M arrays
        !!  of function (derivative) values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: attempt_step => rkv_attempt_step
        !> @brief Perform necessary actions on completion of a successful step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine on_successful_step( &
        !!  class(rk_variable_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The value of the independent variable at the next step.
        !! @param[in] y An N-element array containing the current solution
        !!  values.
        !! @param[in] yn An N-element array containing the solution values at
        !!  the next step.
        procedure, public :: on_successful_step => rkv_on_successful_step
        !> @brief Sets up the interpolation polynomial.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_up_interpolation( &
        !!  class(rk_variable_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) k(:,:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The rk_variable_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The value of the independent variable at the next step.
        !! @param[in] y An N-element array containing the current solution
        !!  values.
        !! @param[in] yn An N-element array containing the solution values at
        !!  the next step.
        !! @param[in] k An N-by-M matrix containing the intermediate step
        !!  function outputs where M is the number of stages of the integrator.
        procedure(rkv_set_up_interp), public, deferred :: set_up_interpolation
        !> @brief Gets the \f$ alpha \f$ control parameter in the PI controller
        !! \f$ h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
        !! e_{n-1}^{\beta}\f$, where \f$ f \f$ is a safety factor, and \f$ k \f$
        !! is the order of the integration method.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_alpha( &
        !!  class(rk_variable_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk_variable_integrator object.
        !! @return The parameter value.
        procedure, public :: get_alpha => rkv_get_alpha
        !> @brief Sets the \f$ alpha \f$ control parameter in the PI controller
        !! \f$ h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
        !! e_{n-1}^{\beta}\f$, where \f$ f \f$ is a safety factor, and \f$ k \f$
        !! is the order of the integration method.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_alpha( &
        !!  class(rk_variable_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_alpha => rkv_set_alpha
        !> @brief Gets the \f$ beta \f$ control parameter in the PI controller
        !! \f$ h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
        !! e_{n-1}^{\beta}\f$, where \f$ f \f$ is a safety factor, and \f$ k \f$
        !! is the order of the integration method.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_beta( &
        !!  class(rk_variable_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rk_variable_integrator object.
        !! @return The parameter value.
        procedure, public :: get_beta => rkv_get_beta
        !> @brief Sets the \f$ beta \f$ control parameter in the PI controller
        !! \f$ h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
        !! e_{n-1}^{\beta}\f$, where \f$ f \f$ is a safety factor, and \f$ k \f$
        !! is the order of the integration method.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_beta( &
        !!  class(rk_variable_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] The parameter value.
        procedure, public :: set_beta => rkv_set_beta
        !> @brief Computes the next step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function compute_next_step_size( &
        !!  class(rk_variable_integrator) this, &
        !!  real(real64) hn, &
        !!  real(real64) en, &
        !!  real(real64) enm1 &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref rk_variable_integrator object.
        !! @param[in] hn The current step size.
        !! @param[in] en The norm of the error for the current step size.
        !! @param[in] enm1 The norm of the error from the previous step size.
        !!
        !! @return The new step size.
        !!
        !! @par Remarks
        !! The step size estimate makes use of a PI type controller such that
        !! \f$ h_{n+1} = f h_n \left( \frac{1}{e_n} \right)^{1/k} e_n^{\alpha} 
        !! e_{n-1}^{\beta}\f$, where \f$ f \f$ is a safety factor, and \f$ k \f$
        !! is the order of the integration method.
        procedure, public :: compute_next_step_size => rkv_next_step
    end type

    interface
        pure function rkv_get_matrix_parameter(this, i, j) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32), intent(in) :: i, j
            real(real64) :: rst
        end function

        pure function rkv_get_array_parameter(this, i) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32), intent(in) :: i
            real(real64) :: rst
        end function

        pure function rkv_get_boolean_parameter(this) result(rst)
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            logical :: rst
        end function

        pure function rkv_get_integer_parameter(this) result(rst)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        subroutine rkv_action(this)
            import rk_variable_integrator
            class(rk_variable_integrator), intent(inout) :: this
        end subroutine

        subroutine rkv_set_up_interp(this, x, xn, y, yn, k)
            use iso_fortran_env
            import rk_variable_integrator
            class(rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
            real(real64), intent(in), dimension(:,:) :: k
        end subroutine

        module subroutine rkv_alloc_workspace(this, neqn, err)
            class(rk_variable_integrator), intent(inout) :: this
            integer(int32), intent(in) :: neqn
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine rkv_reset(this)
            class(rk_variable_integrator), intent(inout) :: this
        end subroutine

        module subroutine rkv_attempt_step(this, sys, h, x, y, yn, en, xprev, &
            yprev, fprev, err)
            class(rk_variable_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn, en
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine rkv_on_successful_step(this, x, xn, y, yn)
            class(rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
        end subroutine

        pure module function rkv_get_alpha(this) result(rst)
            class(rk_variable_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine rkv_set_alpha(this, x)
            class(rk_variable_integrator), intent(inout) :: this
            real(real64) :: x
        end subroutine

        pure module function rkv_get_beta(this) result(rst)
            class(rk_variable_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine rkv_set_beta(this, x)
            class(rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        module function rkv_next_step(this, hn, en, enm1) result(rst)
            class(rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: hn, en, enm1
            real(real64) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the classical Dormand-Prince 4th/5th order integrator.
    !!
    !! @par Example
    !! The following example illustrates how to use this integrator to solve
    !! the forced Duffing model.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq
    !!     use diffeq_models
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     real(real64), parameter :: xmax = 6.0d1
    !!
    !!     ! Local Variables
    !!     type(dprk45_integrator) :: integrator
    !!     type(ode_container) :: mdl
    !!     real(real64), allocatable :: sol(:,:)
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd1, pd2
    !!     class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Define the model
    !!     mdl%fcn => duffing
    !!
    !!     ! Compute the solution
    !!     sol = integrator%solve(mdl, [0.0d0, xmax], [0.0d0, 0.0d0])
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     y2Axis => plt%get_y2_axis()
    !!     lgnd => plt%get_legend()
    !!     call xAxis%set_title("x")
    !!     call yAxis%set_title("y(x)")
    !!     call y2Axis%set_title("y'(x)")
    !!     call plt%set_use_y2_axis(.true.)
    !!     call lgnd%set_is_visible(.true.)
    !!     call xAxis%set_autoscale(.false.)
    !!     call xAxis%set_limits(0.0d0, xmax)
    !!
    !!     call pd1%define_data(sol(:,1), sol(:,2))
    !!     call pd1%set_name("y(x)")
    !!     call plt%push(pd1)
    !!
    !!     call pd2%define_data(sol(:,1), sol(:,3))
    !!     call pd2%set_draw_against_y2(.true.)
    !!     call pd2%set_name("y'(x)")
    !!     call plt%push(pd2)
    !!
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The ODE routine was stored in a seperate module; however, here is the
    !! code for the ODE routine.
    !! @code{.f90}
    !! subroutine duffing(x, y, dydx)
    !!     ! Arguments
    !!     real(real64), intent(in) :: x, y(:)
    !!     real(real64), intent(out) :: dydx(:)
    !!
    !!     ! Model Constants
    !!     real(real64), parameter :: alpha = 1.0d0
    !!     real(real64), parameter :: beta = 5.0d0
    !!     real(real64), parameter :: delta = 2.0d-2
    !!     real(real64), parameter :: gamma = 8.0d0
    !!     real(real64), parameter :: w = 0.5d0
    !!
    !!     ! Equations
    !!     dydx(1) = y(2)
    !!     dydx(2) = gamma * cos(w * x) - delta * y(2) - alpha * y(1) - beta * y(1)**3
    !! end subroutine
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html dprk45_diffing_example_1.png
    type, extends(rk_variable_integrator) :: dprk45_integrator
        logical, private :: m_modelDefined = .false.
        real(real64), private, allocatable, dimension(:,:) :: m_dprk45work
    contains
        !> @brief Defines (initializes) the model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_model(class(dprk45_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref dprk45_integrator object.
        procedure, public :: define_model => dprk45_define_model
        !> @brief Determines if the integrator is an FSAL (first same as last)
        !! integrator (e.g. the 4th/5th order Dormand-Prince integrator).
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function is_fsal(class(dprk45_integrator) this)
        !! @endcode
        !! 
        !! @param[in] this The @ref dprk45_integrator object.
        !! @return Returns true if the integrator is an FSAL type; else,
        !!  returns false.
        procedure, public :: is_fsal => dprk45_is_fsal
        !> @brief Gets the number of stages used by the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_stage_count( &
        !!  class(dprk45_integrator) this &
        !! )
        !! @endcode
        !! 
        !! @param[in] this The @ref dprk45_integrator object.
        !! @return The number of stages.
        procedure, public :: get_stage_count => dprk45_get_stage_count
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order(class(dprk45_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref dprk45_integrator object.
        !! @return The order of the integrator.
        procedure, public :: get_order => dprk45_get_order
        !> @brief Provides interpolation between integration points allowing
        !! for dense output.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine interpolate( &
        !!  class(dprk45_integrator) this, &
        !!  real(real64) xprev, &
        !!  real(real64) xnew, &
        !!  real(real64) x, &
        !!  real(real64) y(:),
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The dprk45_integrator object.
        !! @param[in] xprev The previous value of the independent variable.
        !! @param[in] xnew The updated value of the independent variable.
        !! @param[in] x The value at which to perform the interpolation.
        !! @param[out] y An N-element array containing the interpolated 
        !!  values for each equation.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: interpolate => dprk45_interp
        !> @brief Sets up the interpolation polynomial.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_up_interpolation( &
        !!  class(dprk45_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) k(:,:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The dprk45_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The value of the independent variable at the next step.
        !! @param[in] y An N-element array containing the current solution
        !!  values.
        !! @param[in] yn An N-element array containing the solution values at
        !!  the next step.
        !! @param[in] k An N-by-M matrix containing the intermediate step
        !!  function outputs where M is the number of stages of the integrator.
        procedure, public :: set_up_interpolation => dprk45_set_up_interp
    end type

    ! diffeq_dprk45.f90
    interface
        module subroutine dprk45_define_model(this)
            class(dprk45_integrator), intent(inout) :: this
        end subroutine

        pure module function dprk45_is_fsal(this) result(rst)
            class(dprk45_integrator), intent(in) :: this
            logical :: rst
        end function

        pure module function dprk45_get_stage_count(this) result(rst)
            class(dprk45_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function dprk45_get_order(this) result(rst)
            class(dprk45_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine dprk45_interp(this, xprev, xnew, x, y, err)
            class(dprk45_integrator), intent(in) :: this
            real(real64), intent(in) :: xprev, xnew, x
            real(real64), intent(out), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine dprk45_set_up_interp(this, x, xn, y, yn, k)
            class(dprk45_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
            real(real64), intent(in), dimension(:,:) :: k
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the Bogacki-Shampine 3rd/2nd order integrator.
    !!
    !! @par Example
    !! The following example illustrates how to use this integrator to solve the
    !! Van der Pol model.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use diffeq
    !!     use diffeq_models
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     real(real64), parameter :: xmax = 2.5d1
    !!
    !!     ! Local Variables
    !!     type(bsrk32_integrator) :: integrator
    !!     type(ode_container) :: mdl
    !!     real(real64), allocatable :: sol(:,:)
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: pd
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!
    !!     ! Define the model
    !!     mdl%fcn => vanderpol
    !!
    !!     ! Compute the solution
    !!     sol = integrator%solve(mdl, [0.0d0, xmax], [2.0d0, 0.0d0])
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     call xAxis%set_title("x")
    !!     call yAxis%set_title("y(x)")
    !!     call xAxis%set_autoscale(.false.)
    !!     call xAxis%set_limits(0.0d0, xmax)
    !!
    !!     call pd%define_data(sol(:,1), sol(:,2))
    !!     call pd%set_line_width(2.0)
    !!     call plt%push(pd)
    !!
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The ODE routine was stored in a seperate module; however, here is the
    !! code for the ODE routine.
    !! @code{.f90}
    !! subroutine vanderpol(x, y, dydx)
    !!     ! Arguments
    !!     real(real64), intent(in) :: x, y(:)
    !!     real(real64), intent(out) :: dydx(:)
    !!
    !!     ! Model Constants
    !!     real(real64), parameter :: mu = 5.0d0
    !!
    !!     ! Equations
    !!     dydx(1) = y(2)
    !!     dydx(2) = mu * (1.0d0 - y(1)**2) * y(2) - y(1)
    !! end subroutine
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html bsrk32_vanderpol_example_1.png
    type, extends(rk_variable_integrator) :: bsrk32_integrator
        logical, private :: m_modelDefined = .false.
        real(real64), private, allocatable, dimension(:,:) :: m_bsrk23work
    contains
        !> @brief Defines (initializes) the model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_model(class(bsrk32_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref bsrk32_integrator object.
        procedure, public :: define_model => bsrk32_define_model
        !> @brief Determines if the integrator is an FSAL (first same as last)
        !! integrator (e.g. the 4th/5th order Dormand-Prince integrator).
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function is_fsal(class(bsrk32_integrator) this)
        !! @endcode
        !! 
        !! @param[in] this The @ref bsrk32_integrator object.
        !! @return Returns true if the integrator is an FSAL type; else,
        !!  returns false.
        procedure, public :: is_fsal => bsrk32_is_fsal
        !> @brief Gets the number of stages used by the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_stage_count( &
        !!  class(bsrk32_integrator) this &
        !! )
        !! @endcode
        !! 
        !! @param[in] this The @ref bsrk32_integrator object.
        !! @return The number of stages.
        procedure, public :: get_stage_count => bsrk32_get_stage_count
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order(class(bsrk32_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref bsrk32_integrator object.
        !! @return The order of the integrator.
        procedure, public :: get_order => bsrk32_get_order
        !> @brief Provides interpolation between integration points allowing
        !! for dense output.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine interpolate( &
        !!  class(bsrk32_integrator) this, &
        !!  real(real64) xprev, &
        !!  real(real64) xnew, &
        !!  real(real64) x, &
        !!  real(real64) y(:),
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The bsrk32_integrator object.
        !! @param[in] xprev The previous value of the independent variable.
        !! @param[in] xnew The updated value of the independent variable.
        !! @param[in] x The value at which to perform the interpolation.
        !! @param[out] y An N-element array containing the interpolated 
        !!  values for each equation.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: interpolate => bsrk32_interp
        !> @brief Sets up the interpolation polynomial.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_up_interpolation( &
        !!  class(bsrk32_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) k(:,:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The bsrk32_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The value of the independent variable at the next step.
        !! @param[in] y An N-element array containing the current solution
        !!  values.
        !! @param[in] yn An N-element array containing the solution values at
        !!  the next step.
        !! @param[in] k An N-by-M matrix containing the intermediate step
        !!  function outputs where M is the number of stages of the integrator.
        procedure, public :: set_up_interpolation => bsrk32_set_up_interp
    end type

    ! diffeq_bsrk32.f90
    interface
        module subroutine bsrk32_define_model(this)
            class(bsrk32_integrator), intent(inout) :: this
        end subroutine

        pure module function bsrk32_is_fsal(this) result(rst)
            class(bsrk32_integrator), intent(in) :: this
            logical :: rst
        end function

        pure module function bsrk32_get_stage_count(this) result(rst)
            class(bsrk32_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function bsrk32_get_order(this) result(rst)
            class(bsrk32_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine bsrk32_set_up_interp(this, x, xn, y, yn, k)
            class(bsrk32_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
            real(real64), intent(in), dimension(:,:) :: k
        end subroutine

        module subroutine bsrk32_interp(this, xprev, xnew, x, y, err)
            class(bsrk32_integrator), intent(in) :: this
            real(real64), intent(in) :: xprev, xnew, x
            real(real64), intent(out), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines an implicit Runge-Kutta variable-step integrator.
    type, abstract, extends(rk_variable_integrator) :: &
    implicit_rk_variable_integrator
        ! The most recent successful step size
        real(real64), private :: m_hold = 0.0d0
        ! Is the Jacobian current?
        logical, private :: m_isJacCurrent = .false.
        ! Use a PI step size (true) or a Gustafsson controller (false)
        logical, private :: m_usePI = .true.
    contains
        !> @brief Gets the most recent successful step size.
        !!
        !! @par Synatx
        !! @code{.f90}
        !! real(real64) pure function get_previous_step_size( &
        !!  class(implicit_rk_variable_integrator) this &
        !! )
        !! @endcode
        !! 
        !! @param[in] this The @ref implicit_rk_variable_integrator object.
        !! @return The step size.
        procedure, public :: get_previous_step_size => irk_get_old_step
        !> @brief Sets the most recent successful step size.
        !!
        !! @par Synatx
        !! @code{.f90}
        !! subroutine set_previous_step_size( &
        !!  class(implicit_rk_variable_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !! 
        !! @param[in,out] this The @ref implicit_rk_variable_integrator object.
        !! @param[in] x The step size.
        procedure, public :: set_previous_step_size => irk_set_old_step
        !> @brief Computes the next step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function compute_next_step_size( &
        !!  class(implicit_rk_variable_integrator) this, &
        !!  real(real64) hn, &
        !!  real(real64) en, &
        !!  real(real64) enm1 &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref implicit_rk_variable_integrator object.
        !! @param[in] hn The current step size.
        !! @param[in] en The norm of the error for the current step size.
        !! @param[in] enm1 The norm of the error from the previous step size.
        !!
        !! @return The new step size.
        procedure, public :: compute_next_step_size => &
            irk_compute_next_step_size
        !> @brief Gets a value determining if the Jacobian matrix estimate is
        !! current such that it does not need to be recomputed at this time.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function get_is_jacobian_current( &
        !!  class(implicit_rk_variable_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref implicit_rk_variable_integrator object.
        !! @return True if the Jacobian matrix is current; else, false.
        procedure, public :: get_is_jacobian_current => irk_get_is_jac_current
        !> @brief Sets a value determining if the Jacobian matrix estimate is
        !! current such that it does not need to be recomputed at this time.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_is_jacobian_current( &
        !!  class(implicit_rk_variable_integrator) this, &
        !!  logical x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref implicit_rk_variable_integrator object.
        !! @param[in] x True if the Jacobian matrix is current; else, false.
        procedure, public :: set_is_jacobian_current => irk_set_is_jac_current
        !> @brief Builds the matrix of the form \f$ X = f I - J \f$,
        !! or \f$ X = f M - J \f$ if a mass matrix is defined, and then computes
        !! its LU factorization.  The Jacobian and mass matrices are evaluated
        !! as part of this process, if necessary.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine build_factored_matrix( &
        !!  class(implicit_rk_variable_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref implicit_rk_variable_integrator object.
        !! @param[in,out] sys The @ref ode_container object containing the
        !!  equations to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure(build_factored_newton_matrix_routine), public, deferred :: &
            build_factored_newton_matrix
        !> @brief Gets a parameter determining if a PI step size controller
        !! or a Gustafsson step size controller should be used. The default is
        !! to use a PI step size controller.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function get_use_pi_controller(
        !!  class(implicit_rk_variable_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref implicit_rk_variable_integrator object.
        !! @return True to use a PI controller; else, false to use a
        !!  Gustafsson controller.
        procedure, public :: get_use_pi_controller => irk_get_use_pi_controller
        !> @brief Sets a parameter determining if a PI step size controller
        !! or a Gustafsson step size controller should be used. The default is
        !! to use a PI step size controller.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_use_pi_controller(
        !!  class(implicit_rk_variable_integrator) this, &
        !!  logical x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref implicit_rk_variable_integrator object.
        !! @param[in] x True to use a PI controller; else, false to use a
        !!  Gustafsson controller.
        procedure, public :: set_use_pi_controller => irk_set_use_pi_controller

        ! TO DO: Compute initial step size
    end type

    ! diffeq_implicit_rk.f90
    interface
        subroutine build_factored_newton_matrix_routine(this, sys, h, x, y, err)
            use iso_fortran_env, only : real64
            use ferror, only : errors
            import implicit_rk_variable_integrator
            import ode_container
            class(implicit_rk_variable_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function irk_get_old_step(this) result(rst)
            class(implicit_rk_variable_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine irk_set_old_step(this, x)
            class(implicit_rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        module function irk_compute_next_step_size(this, hn, en, enm1) &
            result(rst)
            class(implicit_rk_variable_integrator), intent(inout) :: this
            real(real64), intent(in) :: hn, en, enm1
            real(real64) :: rst
        end function

        pure module function irk_get_is_jac_current(this) result(rst)
            class(implicit_rk_variable_integrator), intent(in) :: this
            logical :: rst
        end function

        module subroutine irk_set_is_jac_current(this, x)
            class(implicit_rk_variable_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function irk_get_use_pi_controller(this) result(rst)
            class(implicit_rk_variable_integrator), intent(in) :: this
            logical :: rst
        end function

        module subroutine irk_set_use_pi_controller(this, x)
            class(implicit_rk_variable_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a base structure for diagonally implicit Runge-Kutta 
    !! integrators.
    type, abstract, extends(implicit_rk_variable_integrator) :: dirk_integrator
        ! Jacobian matrix workspace
        real(real64), private, allocatable, dimension(:,:) :: m_jac
        ! Mass matrix workspace
        real(real64), private, allocatable, dimension(:,:) :: m_mass
        ! System matrix workspace
        real(real64), private, allocatable, dimension(:,:) :: m_mtx
        ! LU pivot tracking workspace
        integer(int32), private, allocatable, dimension(:) :: m_pvt
        ! NEQN-element workspace array
        real(real64), private, allocatable, dimension(:) :: m_w
        ! NEQN-element workspace array
        real(real64), private, allocatable, dimension(:) :: m_dy
        ! Allowable number of Newton iterations
        integer(int32), private :: m_maxNewtonIter = 7
        ! Newton iteration tolerance
        real(real64), private :: m_newtontol = 1.0d-4
    contains
        !> @brief Initializes the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine initialize( &
        !!  class(dirk_integrator) this, &
        !!  integer(int32) neqn, &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref dirk_integrator object.
        !! @param[in] neqn The number of equations being integrated.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MEMORY_ALLOCATION_ERROR: Occurs if there is a memory 
        !!      allocation issue.
        procedure, public :: initialize => dirk_alloc_workspace
        !> @brief Builds the system matrix of the form \f$ X = f I - J \f$,
        !! or \f$ X = f M - J \f$ if a mass matrix is defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine build_matrix( &
        !!  class(dirk_integrator) this, &
        !!  real(real64) h, &
        !!  real(real64) jac(:,:), &
        !!  real(real64) x(:,:), &
        !!  optional real(real64) m(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref dirk_integrator object.
        !! @param[in] h The current step size.
        !! @param[in] jac The current NEQN-by-NEQN Jacobian matrix.
        !! @param[out] x An NEQN-by-NEQN matrix where the output will be 
        !!  written.
        !! @param[in] m An optional NEQN-by-NEQN mass matrix.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.  Possible errors and warning messages that may be 
        !!  encountered are as follows.
        !!  - DIFFEQ_MATRIX_SIZE_ERROR: Occurs if any of the matrices are not
        !!      sized correctly.
        procedure, public :: build_newton_matrix => dirk_build_matrix
        !> @brief Builds the matrix of the form \f$ X = f I - J \f$,
        !! or \f$ X = f M - J \f$ if a mass matrix is defined, and then computes
        !! its LU factorization.  The Jacobian and mass matrices are evaluated
        !! as part of this process, if necessary.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine build_factored_matrix( &
        !!  class(dirk_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref dirk_integrator object.
        !! @param[in,out] sys The @ref ode_container object containing the
        !!  equations to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: build_factored_newton_matrix => &
            dirk_build_factored_matrix
        !> @brief Gets the maximum allowed number of Newton iterations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_max_newton_iteration_count( &
        !!  class(dirk_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref dirk_integrator object.
        !! @return The iteration limit.
        procedure, public :: get_max_newton_iteration_count => &
            dirk_get_max_newton_iter
        !> @brief Sets the maximum allowed number of Newton iterations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_newton_iteration_count( &
        !!  class(dirk_integrator) this, &
        !!  integer(int32) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref dirk_integrator object.
        !! @param[in] x The iteration limit.
        procedure, public :: set_max_newton_iteration_count => &
            dirk_set_max_newton_iter
        !> @brief Gets the tolerance used to check for convergence of the
        !! Newton iterations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function get_newton_tolerance( &
        !!  class(dirk_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref dirk_integrator object.
        !! @return The tolerance.
        procedure, public :: get_newton_tolerance => dirk_get_newton_tol
        !> @brief Sets the tolerance used to check for convergence of the
        !! Newton iterations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_newton_tolerance( &
        !!  class(dirk_integrator) this, &
        !!  real(real64) x &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref dirk_integrator object.
        !! @param[in] x The tolerance.
        procedure, public :: set_newton_tolerance => dirk_set_newton_tol
        !> @brief Attempts a single integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine attempt_step( &
        !!  class(dirk_integrator) this, &
        !!  class(ode_container) sys, &
        !!  real(real64) h, &
        !!  real(real64) x, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) en(:), &
        !!  optional real(real64) xprev(:), &
        !!  optional real(real64) yprev(:,:), &
        !!  optional real(real64) fprev(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref dirk_integrator object.
        !! @param[in] sys The @ref ode_container object containing the ODEs
        !!  to integrate.
        !! @param[in] h The current step size.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] y An N-element array containing the current values of
        !!  the dependent variables.
        !! @param[out] yn An N-element array where the values of the dependent
        !!  variables at @p x + @p h will be written.
        !! @param[out] en An N-element array where the values of the error
        !!  estimates will be written.
        !! @param[in] xprev An M-element array containing the previous M values
        !!  of the independent variable where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] yprev An M-by-NEQN array containing the previous M arrays
        !!  of dependent variable values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in] fprev An M-by-NEQN array containing the previous M arrays
        !!  of function (derivative) values where M is the order of the method.
        !!  This is typically useful for multi-step methods.  In single-step
        !!  methods this parameter is not used.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: attempt_step => dirk_attempt_step
    end type

    ! diffeq_dirk.f90
    interface
        module subroutine dirk_alloc_workspace(this, neqn, err)
            class(dirk_integrator), intent(inout) :: this
            integer(int32), intent(in) :: neqn
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine dirk_build_matrix(this, h, jac, x, m, err)
            class(dirk_integrator), intent(in) :: this
            real(real64), intent(in) :: h
            real(real64), intent(in), dimension(:,:) :: jac
            real(real64), intent(out), dimension(:,:) :: x
            real(real64), intent(in), dimension(:,:), optional :: m
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine dirk_build_factored_matrix(this, sys, h, x, y, err)
            class(dirk_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine

        pure module function dirk_get_max_newton_iter(this) result(rst)
            class(dirk_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine dirk_set_max_newton_iter(this, x)
            class(dirk_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine

        pure module function dirk_get_newton_tol(this) result(rst)
            class(dirk_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine dirk_set_newton_tol(this, x)
            class(dirk_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        module subroutine dirk_attempt_step(this, sys, h, x, y, yn, en, &
            xprev, yprev, fprev, err)
            class(dirk_integrator), intent(inout) :: this
            class(ode_container), intent(inout) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn, en
            real(real64), intent(in), optional, dimension(:) :: xprev
            real(real64), intent(in), optional, dimension(:,:) :: yprev
            real(real64), intent(inout), optional, dimension(:,:) :: fprev
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines a diagonally implicit 4th order Runge-Kutta integrator
    !! suitable for integrating stiff systems of differential equations.
    type, extends(dirk_integrator) :: sdirk4_integrator
        logical, private :: m_modelDefined = .false.
    contains
        !> @brief Returns the order of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_order( &
        !!  class(sdirk4_integrator) this &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref sdirk4_integrator object.
        !! @return The order of the integrator.
        procedure, public :: get_order => sd4_get_order
        !> @brief Gets the number of stages used by the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! integer(int32) pure function get_stage_count( &
        !!  class(sdirk4_integrator) this &
        !! )
        !! @endcode
        !! 
        !! @param[in] this The @ref sdirk4_integrator object.
        !! @return The number of stages.
        procedure, public :: get_stage_count => sd4_get_stage_count
        !> @brief Defines (initializes) the model parameters.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_model(class(sdirk4_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The @ref sdirk4_integrator object.
        procedure, public :: define_model => sd4_define_model
        !> @brief Determines if the integrator is an FSAL (first same as last)
        !! integrator (e.g. the 4th/5th order Dormand-Prince integrator).
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical pure function is_fsal(class(sdirk4_integrator) this)
        !! @endcode
        !! 
        !! @param[in] this The @ref sdirk4_integrator object.
        !! @return Returns true if the integrator is an FSAL type; else,
        !!  returns false.
        procedure, public :: is_fsal => sd4_is_fsal
        !> @brief Provides interpolation between integration points allowing
        !! for dense output.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine interpolate( &
        !!  class(sdirk4_integrator) this, &
        !!  real(real64) xprev, &
        !!  real(real64) xnew, &
        !!  real(real64) x, &
        !!  real(real64) y(:),
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in] this The sdirk4_integrator object.
        !! @param[in] xprev The previous value of the independent variable.
        !! @param[in] xnew The updated value of the independent variable.
        !! @param[in] x The value at which to perform the interpolation.
        !! @param[out] y An N-element array containing the interpolated 
        !!  values for each equation.
        !! @param[in,out] err An optional errors-based object that if provided 
        !!  can be used to retrieve information relating to any errors 
        !!  encountered during execution. If not provided, a default 
        !!  implementation of the errors class is used internally to provide 
        !!  error handling.
        procedure, public :: interpolate => sd4_interp
        !> @brief Sets up the interpolation polynomial.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_up_interpolation( &
        !!  class(sdirk4_integrator) this, &
        !!  real(real64) x, &
        !!  real(real64) xn, &
        !!  real(real64) y(:), &
        !!  real(real64) yn(:), &
        !!  real(real64) k(:,:) &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The sdirk4_integrator object.
        !! @param[in] x The current value of the independent variable.
        !! @param[in] xn The value of the independent variable at the next step.
        !! @param[in] y An N-element array containing the current solution
        !!  values.
        !! @param[in] yn An N-element array containing the solution values at
        !!  the next step.
        !! @param[in] k An N-by-M matrix containing the intermediate step
        !!  function outputs where M is the number of stages of the integrator.
        procedure, public :: set_up_interpolation => sd4_set_up_interp
    end type

    ! diffeq_sdirk4.f90
    interface
        pure module function sd4_get_order(this) result(rst)
            class(sdirk4_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        pure module function sd4_get_stage_count(this) result(rst)
            class(sdirk4_integrator), intent(in) :: this
            integer(int32) :: rst
        end function

        module subroutine sd4_define_model(this)
            class(sdirk4_integrator), intent(inout) :: this
        end subroutine

        pure module function sd4_is_fsal(this) result(rst)
            class(sdirk4_integrator), intent(in) :: this
            logical :: rst
        end function

        module subroutine sd4_set_up_interp(this, x, xn, y, yn, k)
            class(sdirk4_integrator), intent(inout) :: this
            real(real64), intent(in) :: x, xn
            real(real64), intent(in), dimension(:) :: y, yn
            real(real64), intent(in), dimension(:,:) :: k
        end subroutine

        module subroutine sd4_interp(this, xprev, xnew, x, y, err)
            class(sdirk4_integrator), intent(in) :: this
            real(real64), intent(in) :: xprev, xnew, x
            real(real64), intent(out), dimension(:) :: y
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
end module