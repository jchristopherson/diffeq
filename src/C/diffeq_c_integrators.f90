! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module diffeq_c_integrators
    use iso_fortran_env
    use iso_c_binding
    use diffeq_c_containers
    use diffeq_c_errors
    use diffeq
    implicit none
    private
    public :: c_ode
    public :: c_ode_jacobian
    public :: c_ode_mass_matrix
    public :: c_ode_settings
    public :: c_default_ode_settings
    public :: c_runge_kutta_23
    public :: c_runge_kutta_45
    public :: c_runge_kutta_853
    public :: c_tsitouras
    public :: c_adams
    public :: c_rosenbrock
    public :: c_kennedy_carpenter_4
    public :: c_kennedy_carpenter_5
    public :: c_bdf

    integer, parameter :: C_RK23 = 1
    integer, parameter :: C_RK45 = 2
    integer, parameter :: C_RK853 = 3
    integer, parameter :: C_TSIT = 4
    integer, parameter :: C_RBRK = 5
    integer, parameter :: C_KCARP4 = 6
    integer, parameter :: C_KCARP5 = 7
    integer, parameter :: C_BACKDIFF = 8
    integer, parameter :: C_ADMS = 9

    interface
        subroutine c_ode(n, x, y, dydx) bind(C)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: dydx(n)
        end subroutine

        subroutine c_ode_jacobian(m, n, x, y, jac, ldj) bind(C)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: m
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: jac(ldj,n)
            integer(c_int), intent(in), value :: ldj
        end subroutine

        subroutine c_ode_mass_matrix(n, x, y, m, ldm) bind(C)
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
            real(c_double), intent(in), value :: x
            real(c_double), intent(in) :: y(n)
            real(c_double), intent(out) :: m(ldm,n)
            integer(c_int), intent(in), value :: ldm
        end subroutine
    end interface

    type, bind(C) :: c_ode_settings
        !! A container for ODE settings.
        logical(c_bool) :: mass_matrix_constant
            !! True if the mass matrix is constant regardless of state, or
            !! false if the mass matrix is state-dependent.
        real(c_double) :: minimum_step_size
            !! The minimum allowed integrator step size.
        real(c_double) :: maximum_step_size
            !! The maximum allowed integrator step size.
        integer(c_int) :: step_limit
            !! The limit on the maximum total number of integration steps.
        logical(c_bool) :: allow_overshoot
            !! True if the solver is allowed to overshoot the final value in
            !! the integration range and interpolate back; else, false if the
            !! solver must terminate on the final value.
        real(c_double) :: beta
            !! PI step size controller exponent.  Set to zero to turn off PI
            !! step size control.
        real(c_double) :: absolute_tolerance
            !! The absolute tolerance applied to each equation.
        real(c_double) :: relative_tolerance
            !! The relative tolerance applied to each equation.
        integer(c_int) :: newton_step_limit
            !! The maximum number of Newton iterations allowed per implicit
            !! solve (BDF and Kennedy-Carpenter).
        real(c_double) :: newton_tolerance
            !! The convergence tolerance for Newton iteration (BDF and
            !! Kennedy-Carpenter).
    end type

contains
! ------------------------------------------------------------------------------
    subroutine c_default_ode_settings(settings) &
        bind(C, name = "c_default_ode_settings")
        !! Fills the supplied [[c_ode_settings]] object with default settings.
        type(c_ode_settings), intent(inout) :: settings
            !! The [[c_ode_settings]] object.
        settings%absolute_tolerance = 1.0d-6
        settings%allow_overshoot = .true.
        settings%beta = 0.0d0
        settings%mass_matrix_constant = .true.
        settings%maximum_step_size = huge(1.0d0)
        settings%minimum_step_size = 1.0d1 * epsilon(1.0d0)
        settings%relative_tolerance = 1.0d-6
        settings%step_limit = 1000000
        settings%newton_step_limit = 12
        settings%newton_tolerance = 1.0d-8
    end subroutine

! ------------------------------------------------------------------------------
    pure function c_valid_newton_settings(settings) result(rst)
        !! Determines if the Newton iteration settings are usable.
        type(c_ode_settings), intent(in) :: settings
            !! The [[c_ode_settings]] object.
        logical :: rst
            !! True if the settings are valid; else, false.
        rst = settings%newton_step_limit >= 1 .and. &
            settings%newton_tolerance > 0.0d0 .and. &
            settings%newton_tolerance <= huge(1.0d0)
    end function

! ------------------------------------------------------------------------------
    function c_integrate_ode(method, nx, nstates, fcn, x, initcond, &
        settings, output, jac, mass) result(rst)
        integer(c_int), intent(in) :: method
            !! The integrator to use.
        integer(c_int), intent(in) :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in) :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        procedure(c_ode), intent(in), pointer :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        procedure(c_ode_jacobian), intent(in), pointer, optional :: jac
            !! A pointer to the Jacobian calculation routine.  If not supplied,
            !! a finite difference approach is utilized if a Jacobian matrix
            !! is needed.
        procedure(c_ode_mass_matrix), intent(in), pointer, optional :: mass
            !! A pointer to the mass matrix routine.  If not supplied, no 
            !! mass matrix is utilized.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        class(ode_integrator), pointer :: integrator
        type(runge_kutta_23), target :: rk23
        type(runge_kutta_45), target :: rk45
        type(runge_kutta_853), target :: rk853
        type(tsitouras_54), target :: tsit
        type(rosenbrock), target :: rsbrk
        type(kennedy_carpenter_4), target :: kc4
        type(kennedy_carpenter_5), target :: kc5
        type(bdf), target :: bdf_i
        type(adams), target :: adams_i
        type(ode_container) :: model
        logical :: mass_matrix_dependent
        type(ode_results) :: ode_output

        rst = C_NO_ERROR
        select case (method)
        case (C_RK23)
            integrator => rk23
        case (C_RK45)
            integrator => rk45
        case (C_RK853)
            integrator => rk853
        case (C_TSIT)
            integrator => tsit
        case (C_RBRK)
            integrator => rsbrk
        case (C_KCARP4)
            if (.not.c_valid_newton_settings(settings)) then
                rst = C_INVALID_INPUT_ERROR
                return
            end if
            call kc4%set_newton_step_limit(settings%newton_step_limit)
            call kc4%set_newton_tolerance(settings%newton_tolerance)
            integrator => kc4
        case (C_KCARP5)
            if (.not.c_valid_newton_settings(settings)) then
                rst = C_INVALID_INPUT_ERROR
                return
            end if
            call kc5%set_newton_step_limit(settings%newton_step_limit)
            call kc5%set_newton_tolerance(settings%newton_tolerance)
            integrator => kc5
        case (C_BACKDIFF)
            call bdf_i%set_newton_step_limit(settings%newton_step_limit)
            call bdf_i%set_newton_tolerance(settings%newton_tolerance)
            integrator => bdf_i
        case (C_ADMS)
            integrator => adams_i
        case default
            rst = C_INVALID_INPUT_ERROR
            return
        end select

        if (nx < 2) then
            rst = C_INVALID_INPUT_ERROR
            return
        end if
        if (nstates < 1) then
            rst = C_INVALID_INPUT_ERROR
            return
        end if

        if (present(jac)) model%jacobian => internal_jac
        if (present(mass)) model%mass_matrix => internal_mass
        model%fcn => internal_ode
        mass_matrix_dependent = .not.settings%mass_matrix_constant
        call model%set_is_mass_matrix_dependent(mass_matrix_dependent)

        call integrator%set_minimum_step_size(settings%minimum_step_size)
        call integrator%set_maximum_step_size(settings%maximum_step_size)
        call integrator%set_step_limit(settings%step_limit)
        call integrator%set_allow_overshoot(logical(settings%allow_overshoot))
        call integrator%set_absolute_tolerance(settings%absolute_tolerance)
        call integrator%set_relative_tolerance(settings%relative_tolerance)
        call integrator%set_step_size_control_parameter(settings%beta)

        call integrator%solve(model, x, initcond)
        ode_output = integrator%get_ode_results()
        call c_create_ode_results(output, ode_output)

    contains
        subroutine internal_ode(x_, y_, dydx_, args_)
            real(real64), intent(in) :: x_
            real(real64), intent(in), dimension(:) :: y_
            real(real64), intent(out), dimension(:) :: dydx_
            class(*), intent(inout), optional :: args_
            call fcn(size(y_), x_, y_, dydx_)
        end subroutine

        subroutine internal_jac(x_, y_, jac_, args_)
            real(real64), intent(in) :: x_
            real(real64), intent(in), dimension(:) :: y_
            real(real64), intent(out), dimension(:,:) :: jac_
            class(*), intent(inout), optional :: args_
            call jac(size(jac_, 1), size(jac_, 2), x_, y_, jac_, size(jac_, 1))
        end subroutine

        subroutine internal_mass(x_, y_, m_, args_)
            real(real64), intent(in) :: x_
            real(real64), intent(in), dimension(:) :: y_
            real(real64), intent(out), dimension(:,:) :: m_
            class(*), intent(inout), optional :: args_
            call mass(size(m_, 1), x_, y_, m_, size(m_, 1))
        end subroutine
    end function

! ------------------------------------------------------------------------------
    function c_runge_kutta_23(nx, nstates, fcn, x, initcond, settings, output) &
        result(rst) bind(C, name = "c_runge_kutta_23")
        !! Utilizes the 3rd order Bogacki-Shampine Runge-Kutta integrator to 
        !! solve a system of first order differential equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)
        rst = c_integrate_ode(C_RK23, nx, nstates, fcn_ptr, x, initcond, &
            settings, output)
    end function

! ------------------------------------------------------------------------------
    function c_runge_kutta_45(nx, nstates, fcn, x, initcond, settings, output) &
        result(rst) bind(C, name = "c_runge_kutta_45")
        !! Utilizes the 5th order Dormand-Prince Runge-Kutta integrator to 
        !! solve a system of first order differential equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)
        rst = c_integrate_ode(C_RK45, nx, nstates, fcn_ptr, x, initcond, &
            settings, output)
    end function

! ------------------------------------------------------------------------------
    function c_runge_kutta_853(nx, nstates, fcn, x, initcond, settings, output) &
        result(rst) bind(C, name = "c_runge_kutta_853")
        !! Utilizes the 8th order Hairer Runge-Kutta integrator to 
        !! solve a system of first order differential equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)
        rst = c_integrate_ode(C_RK853, nx, nstates, fcn_ptr, x, initcond, &
            settings, output)
    end function

! ------------------------------------------------------------------------------
    function c_tsitouras(nx, nstates, fcn, x, initcond, settings, output) &
        result(rst) bind(C, name = "c_tsitouras")
        !! Utilizes the 5th order Tsitouras Runge-Kutta integrator to 
        !! solve a system of first order differential equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)
        rst = c_integrate_ode(C_TSIT, nx, nstates, fcn_ptr, x, initcond, &
            settings, output)
    end function

! ------------------------------------------------------------------------------
    function c_adams(nx, nstates, fcn, x, initcond, settings, mass, output) &
        result(rst) bind(C, name = "c_adams")
        !! Utilizes a variable order Adams predictor-corrector integrator to 
        !! solve a system of first order differential equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_funptr), value, intent(in) :: mass
            !! A pointer to the routine for computing the mass matrix.  Pass
            !! NULL if no routine is necessary.  This routine cannot handle
            !! singular mass matrices.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        logical :: use_mass
        procedure(c_ode), pointer :: fcn_ptr
        procedure(c_ode_mass_matrix), pointer :: mass_ptr

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)

        use_mass = c_associated(mass)
        if (use_mass) call c_f_procpointer(mass, mass_ptr)

        if (use_mass) then
            rst = c_integrate_ode(C_ADMS, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, mass = mass_ptr)
        else
            rst = c_integrate_ode(C_ADMS, nx, nstates, fcn_ptr, x, initcond, &
                settings, output)
        end if
    end function

! ------------------------------------------------------------------------------
    function c_rosenbrock(nx, nstates, fcn, x, initcond, settings, jac, mass, &
        output) result(rst) bind(C, name = "c_rosenbrock")
        !! Utilizes a 4th order Rosenbrock integrator to solve a system of 
        !! first order differential equations or index-1 differential 
        !! algebraic equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_funptr), value, intent(in) :: jac
            !! A pointer to the routine for computing the Jacobian matrix.
            !! Pass NULL if no routine is necessary and allowing the integrator
            !! to use a finite-difference estimate is acceptable.
        type(c_funptr), value, intent(in) :: mass
            !! A pointer to the routine for computing the mass matrix.  Pass
            !! NULL if no routine is necessary.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr
        procedure(c_ode_jacobian), pointer :: jac_ptr
        procedure(c_ode_mass_matrix), pointer :: mass_ptr
        logical :: use_jac, use_mass

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)

        use_jac = c_associated(jac)
        if (use_jac) call c_f_procpointer(jac, jac_ptr)

        use_mass = c_associated(mass)
        if (use_mass) call c_f_procpointer(mass, mass_ptr)

        if (use_jac .and. use_mass) then
            rst = c_integrate_ode(C_RBRK, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr, mass = mass_ptr)
        else if (use_jac .and. .not.use_mass) then
            rst = c_integrate_ode(C_RBRK, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr)
        else if (.not.use_jac .and. use_mass) then
            rst = c_integrate_ode(C_RBRK, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, mass = mass_ptr)
        else
            rst = c_integrate_ode(C_RBRK, nx, nstates, fcn_ptr, x, initcond, &
                settings, output)
        end if
    end function

! ------------------------------------------------------------------------------
    function c_kennedy_carpenter_4(nx, nstates, fcn, x, initcond, settings, &
        jac, mass, output) result(rst) bind(C, name = "c_kennedy_carpenter_4")
        !! Utilizes a 4th order Kennedy-Carpenter integrator to solve a system 
        !! of first order differential equations or index-1 differential 
        !! algebraic equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_funptr), value, intent(in) :: jac
            !! A pointer to the routine for computing the Jacobian matrix.
            !! Pass NULL if no routine is necessary and allowing the integrator
            !! to use a finite-difference estimate is acceptable.
        type(c_funptr), value, intent(in) :: mass
            !! A pointer to the routine for computing the mass matrix.  Pass
            !! NULL if no routine is necessary.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr
        procedure(c_ode_jacobian), pointer :: jac_ptr
        procedure(c_ode_mass_matrix), pointer :: mass_ptr
        logical :: use_jac, use_mass

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)

        use_jac = c_associated(jac)
        if (use_jac) call c_f_procpointer(jac, jac_ptr)

        use_mass = c_associated(mass)
        if (use_mass) call c_f_procpointer(mass, mass_ptr)

        if (use_jac .and. use_mass) then
            rst = c_integrate_ode(C_KCARP4, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr, mass = mass_ptr)
        else if (use_jac .and. .not.use_mass) then
            rst = c_integrate_ode(C_KCARP4, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr)
        else if (.not.use_jac .and. use_mass) then
            rst = c_integrate_ode(C_KCARP4, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, mass = mass_ptr)
        else
            rst = c_integrate_ode(C_KCARP4, nx, nstates, fcn_ptr, x, initcond, &
                settings, output)
        end if
    end function

! ------------------------------------------------------------------------------
    function c_kennedy_carpenter_5(nx, nstates, fcn, x, initcond, settings, &
        jac, mass, output) result(rst) bind(C, name = "c_kennedy_carpenter_5")
        !! Utilizes a 5th order Kennedy-Carpenter integrator to solve a system 
        !! of first order differential equations or index-1 differential 
        !! algebraic equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_funptr), value, intent(in) :: jac
            !! A pointer to the routine for computing the Jacobian matrix.
            !! Pass NULL if no routine is necessary and allowing the integrator
            !! to use a finite-difference estimate is acceptable.
        type(c_funptr), value, intent(in) :: mass
            !! A pointer to the routine for computing the mass matrix.  Pass
            !! NULL if no routine is necessary.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr
        procedure(c_ode_jacobian), pointer :: jac_ptr
        procedure(c_ode_mass_matrix), pointer :: mass_ptr
        logical :: use_jac, use_mass

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)

        use_jac = c_associated(jac)
        if (use_jac) call c_f_procpointer(jac, jac_ptr)

        use_mass = c_associated(mass)
        if (use_mass) call c_f_procpointer(mass, mass_ptr)

        if (use_jac .and. use_mass) then
            rst = c_integrate_ode(C_KCARP5, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr, mass = mass_ptr)
        else if (use_jac .and. .not.use_mass) then
            rst = c_integrate_ode(C_KCARP5, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, jac = jac_ptr)
        else if (.not.use_jac .and. use_mass) then
            rst = c_integrate_ode(C_KCARP5, nx, nstates, fcn_ptr, x, initcond, &
                settings, output, mass = mass_ptr)
        else
            rst = c_integrate_ode(C_KCARP5, nx, nstates, fcn_ptr, x, initcond, &
                settings, output)
        end if
    end function

! ------------------------------------------------------------------------------
    function c_bdf(nx, nstates, fcn, x, initcond, settings, jac, mass, &
        output) result(rst) bind(C, name = "c_bdf")
        !! Utilizes a variable order BDF integrator to solve a system of 
        !! first order differential equations or index-1 differential 
        !! algebraic equations.
        integer(c_int), intent(in), value :: nx
            !! The number of solution points.  This value must be at least 2;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        integer(c_int), intent(in), value :: nstates
            !! The number of state equations.  This value must be at least 1;
            !! else, a C_INVALID_INPUT_ERROR is returned.
        type(c_funptr), value, intent(in) :: fcn
            !! A pointer to the function containing the equations to integrate.
        real(c_double), intent(in) :: x(nx)
            !! An NX-element array containing the desired solution points.
        real(c_double), intent(in) :: initcond(nstates)
            !! An NSTATES-element array containing the initial conditions.
        type(c_ode_settings), intent(in) :: settings
            !! The integrator settings.
        type(c_funptr), value, intent(in) :: jac
            !! A pointer to the routine for computing the Jacobian matrix.
            !! Pass NULL if no routine is necessary and allowing the integrator
            !! to use a finite-difference estimate is acceptable.
        type(c_funptr), value, intent(in) :: mass
            !! A pointer to the routine for computing the mass matrix.  Pass
            !! NULL if no routine is necessary.
        type(c_ptr), intent(out) :: output
            !! The output c_ode_results object.  Memory must be freed by 
            !! c_free_ode_results.
        integer(c_int) :: rst
            !! An error flag denoting the status of the operation.

        procedure(c_ode), pointer :: fcn_ptr
        procedure(c_ode_jacobian), pointer :: jac_ptr
        procedure(c_ode_mass_matrix), pointer :: mass_ptr
        logical :: use_jac, use_mass

        if (.not.c_associated(fcn)) then
            rst = C_NULL_POINTER_ERROR
            return
        end if
        call c_f_procpointer(fcn, fcn_ptr)

        use_jac = c_associated(jac)
        if (use_jac) call c_f_procpointer(jac, jac_ptr)

        use_mass = c_associated(mass)
        if (use_mass) call c_f_procpointer(mass, mass_ptr)

        if (use_jac .and. use_mass) then
            rst = c_integrate_ode(C_BACKDIFF, nx, nstates, fcn_ptr, x, &
                initcond, settings, output, jac = jac_ptr, mass = mass_ptr)
        else if (use_jac .and. .not.use_mass) then
            rst = c_integrate_ode(C_BACKDIFF, nx, nstates, fcn_ptr, x, &
                initcond, settings, output, jac = jac_ptr)
        else if (.not.use_jac .and. use_mass) then
            rst = c_integrate_ode(C_BACKDIFF, nx, nstates, fcn_ptr, x, &
                initcond, settings, output, mass = mass_ptr)
        else
            rst = c_integrate_ode(C_BACKDIFF, nx, nstates, fcn_ptr, x, &
                initcond, settings, output)
        end if
    end function

! ------------------------------------------------------------------------------
end module