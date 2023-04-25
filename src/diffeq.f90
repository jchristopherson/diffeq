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
    public :: variable_step_integrator
    public :: DIFFEQ_MEMORY_ALLOCATION_ERROR
    public :: DIFFEQ_NULL_POINTER_ERROR
    public :: DIFFEQ_MATRIX_SIZE_ERROR
    public :: DIFFEQ_ARRAY_SIZE_ERROR
    public :: DIFFEQ_INVALID_INPUT_ERROR

! ------------------------------------------------------------------------------
    integer(int32), parameter :: DIFFEQ_MEMORY_ALLOCATION_ERROR = 10000
    integer(int32), parameter :: DIFFEQ_NULL_POINTER_ERROR = 10001
    integer(int32), parameter :: DIFFEQ_MATRIX_SIZE_ERROR = 10002
    integer(int32), parameter :: DIFFEQ_ARRAY_SIZE_ERROR = 10003
    integer(int32), parameter :: DIFFEQ_INVALID_INPUT_ERROR = 10004

! ------------------------------------------------------------------------------
    !> @brief A container for the routine containing the ODEs to integrate.
    type ode_container
    private
        ! A value determining if the mass matrix is state dependent such that it
        ! must be recomputed at each step.
        logical :: m_massDependent = .true.
        ! Jacobian calculation workspace array.
        real(real64), allocatable, dimension(:) :: m_jwork
        ! Finite difference step size.
        real(real64) :: m_fdStep = sqrt(epsilon(1.0d0))
        !> @brief A pointer to the routine containing the ODEs to integrate.
        procedure(ode), pointer, public, nopass :: fcn => null()
        !> @brief A pointer to the routine containing the analytical Jacobian.
        !! If supplied, this routine is utilized; however, if null, a finite
        !! difference approximation is utilized.
        procedure(ode_jacobian), pointer, public, nopass :: jacobian => null()
        !> @brief A pointer to the routine containing the mass matrix for the
        !! system.  If set to null (the default), an identity mass matrix will
        !! be assumed.
        procedure(ode_mass_matrix), pointer, public, nopass :: mass_matrix => null()
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
        !!  real(real64) dydx(:), &
        !!  real(real64) jac(:,:), &
        !!  optional class(errors) err &
        !! )
        !! @endcode
        !!
        !! @param[in,out] this The @ref ode_container object.
        !! @param[in] x The current independent variable value.
        !! @param[in] y An N-element array containing the current dependent
        !!  variable values.
        !! @param[in] dydx An N-element array containing the current derivative
        !!  values.
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
    end type

    ! diffeq_ode_container.f90
    interface
        subroutine ode(x, y, dydx)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        subroutine ode_jacobian(x, y, dydx, jac)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y, dydx
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        subroutine ode_mass_matrix(x, y, dydx, m)
            use iso_fortran_env
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y, dydx
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

        module subroutine oc_jacobian(this, x, y, dydx, jac, err)
            class(ode_container), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y, dydx
            real(real64), intent(out), dimension(:,:) :: jac
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief The most basic ODE integrator object capable of integrating
    !! systems of ODE's.
    type, abstract :: ode_integrator
    private
    contains
        procedure(ode_solver), public, pass, deferred :: solve
        procedure(ode_integer_inquiry), public, pass, deferred :: get_order
    end type

    interface
        function ode_solver(this, sys, x, iv, err) result(rst)
            use iso_fortran_env
            use ferror
            import ode_integrator
            import ode_container
            class(ode_integrator), intent(inout) :: this
            class(ode_container), intent(in) :: sys
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
        procedure, public :: solve => fsi_solver
        procedure(ode_fixed_step), public, pass, deferred :: step
    end type

    ! diffeq_fs_integrator.f90
    interface
        subroutine ode_fixed_step(this, sys, h, x, y, yn, err)
            use iso_fortran_env 
            use ferror   
            import fixed_step_integrator
            import ode_container
            class(fixed_step_integrator), intent(inout) :: this
            class(ode_container), intent(in) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module function fsi_solver(this, sys, x, iv, err) result(rst)
            class(fixed_step_integrator), intent(inout) :: this
            class(ode_container), intent(in) :: sys
            real(real64), intent(in), dimension(:) :: x, iv
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines an explicit, Runge-Kutta fixed-step integrator.
    type, abstract, extends(fixed_step_integrator) :: rk_fixed_integrator
    private
        ! Workspace matrix
        real(real64), allocatable, dimension(:,:) :: m_work
    contains
        ! Use to allocate internal workspaces.  This routine only takes action
        ! if the workspace array(s) are not sized properly for the application.
        procedure, private :: allocate_workspace => rkf_alloc_workspace
        procedure, public :: step => rkf_step
        procedure(rkf_get_matrix_parameter), deferred, public :: get_a_parameter
        procedure(rkf_get_array_parameter), deferred, public :: get_b_parameter
        procedure(rkf_get_array_parameter), deferred, public :: get_c_parameter
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

        module subroutine rkf_step(this, sys, h, x, y, yn, err)
            class(rk_fixed_integrator), intent(inout) :: this
            class(ode_container), intent(in) :: sys
            real(real64), intent(in) :: h, x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: yn
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    !> @brief Defines the explicit, 4th order, Runge-Kutta fixed-step 
    !! integrator.
    type, extends(rk_fixed_integrator) :: rk4_fixed_integrator
    contains
        procedure, public :: get_order => rk4_get_order
        procedure, public :: get_a_parameter => rk4_get_a_param
        procedure, public :: get_b_parameter => rk4_get_b_param
        procedure, public :: get_c_parameter => rk4_get_c_param
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
    !> @brief Defines a variable-step integrator.
    type, abstract, extends(ode_integrator) :: variable_step_integrator
    private
        real(real64) :: m_safetyfactor = 0.9d0
        real(real64) :: m_alpha = 0.7d0
        real(real64) :: m_beta = 0.4d0
        real(real64) :: m_maxstep = huge(1.0d0)
    contains
        procedure, public :: get_safety_factor => vsi_get_safety_factor
        procedure, public :: set_safety_factor => vsi_set_safety_factor
        procedure, public :: get_alpha => vsi_get_alpha
        procedure, public :: set_alpha => vsi_set_alpha
        procedure, public :: get_beta => vsi_get_beta
        procedure, public :: set_beta => vsi_set_beta
        procedure, public :: get_max_step_size => vsi_get_max_step
        procedure, public :: set_max_step_size => vsi_set_max_step
        !> @brief Computes a normalized estimates the error for the given step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function estimate_error( &
        !!  class(variable_step_integrator) this, &
        !!  real(real64) y(:), &
        !!  real(real64) ys(:), &
        !!  real(real64) atol(:), &
        !!  real(real64) rtol(:) &
        !! )
        !! @endcode
        !!
        !! @param[in] this The variable_step_integrator object.
        !! @param[in] y The current solution estimate.
        !! @param[in] ys A suplemental solution estimate typically the result
        !!  of an embedded estimate.
        !! @param[in] atol The absolute tolerance value.
        !! @param[in] rtol The relative tolerance value.
        !!
        !! @return The normalized error estimate.
        procedure, public :: estimate_error => vsi_estimate_error
        procedure, public :: compute_next_step_size => vsi_next_step
    end type

    ! diffeq_vs_integrator.f90
    interface
        pure module function vsi_get_safety_factor(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_safety_factor(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function vsi_get_alpha(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_alpha(this, x)
            class(variable_step_integrator), intent(inout) :: this
            real(real64) :: x
        end subroutine

        pure module function vsi_get_beta(this) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64) :: rst
        end function

        module subroutine vsi_set_beta(this, x)
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

        pure module function vsi_estimate_error(this, y, ys, atol, rtol) &
            result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64), intent(in), dimension(:) :: y, ys, atol, rtol
            real(real64) :: rst
        end function

        pure module function vsi_next_step(this, hn, en, enm1) result(rst)
            class(variable_step_integrator), intent(in) :: this
            real(real64), intent(in) :: hn, en, enm1
            real(real64) :: rst
        end function
    end interface
end module