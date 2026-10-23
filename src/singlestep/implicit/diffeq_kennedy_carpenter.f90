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

module diffeq_kennedy_carpenter
    use iso_fortran_env
    use diffeq_base
    use diffeq_singlestep
    use diffeq_errors
    use linalg
    implicit none
    private
    public :: kennedy_carpenter
    public :: kennedy_carpenter_4
    public :: kennedy_carpenter_5

    type, abstract, extends(single_step_integrator) :: kennedy_carpenter
        !! Shared implementation for Kennedy--Carpenter ESDIRK methods.
        !!
        !! The stage equations are solved with Newton iteration.  A supplied
        !! mass matrix is incorporated in each stage system as
        !! \f$M-h a_{ii}J\f$, while the embedded tableau supplies the error
        !! estimate used by the inherited adaptive driver.  These solvers use
        !! the inherited configurable PI step-size controller.
        real(real64), private, allocatable, dimension(:,:) :: m_mass
        real(real64), private, allocatable, dimension(:) :: m_startDerivative
        logical, private :: m_massComputed = .false.
    contains
        procedure, public :: solve => kc_solve
        procedure, public :: pre_step_action => kc_pre_step
            !! Performs any pre-step actions.
        procedure, public :: attempt_step => kc_attempt_step
            !! Attempts an integration step for this integrator.
        procedure, public :: post_step_action => kc_post_step
            !! Performs any post-step actions.
        procedure, public :: interpolate => kc_interpolate
            !! Performs the interpolation.
        procedure, public :: get_is_fsal => kc_get_is_fsal
            !! Gets a logical parameter stating if this is a first-same-as-last
            !! (FSAL) integrator.
        procedure, public :: get_stage_count => kc_get_stage_count
            !! Gets the stage count for this integrator.
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_4
        !! Fourth-order ARK4(3)6L[2]SA ESDIRK method.
        !! The inherited PI step-size controller uses the embedded error
        !! estimate to select the next step size.
    contains
        procedure, public :: get_order => kc4_get_order
            !! Gets the order of the integrator.
    end type

    type, extends(kennedy_carpenter) :: kennedy_carpenter_5
        !! Fifth-order ARK5(4)8L[2]SA ESDIRK method.
        !! The inherited PI step-size controller uses the embedded error
        !! estimate to select the next step size.
    contains
        procedure, public :: get_order => kc5_get_order
            !! Gets the order of the integrator.
    end type

contains
! ------------------------------------------------------------------------------
subroutine kc_solve(this, sys, x, iv, args)
    class(kennedy_carpenter), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: x, iv
    class(*), intent(inout), optional, target :: args

    this%m_massComputed = .false.
    if (allocated(this%m_mass)) deallocate(this%m_mass)
    if (allocated(this%m_startDerivative)) deallocate(this%m_startDerivative)
    call ssi_ode_solver(this, sys, x, iv, args)
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_pre_step(this, prevs, sys, h, x, y, f, args)
    !! Placeholder routine for any pre-step actions.
    !!
    !! The Jacobian and mass matrices are state-dependent within each stage,
    !! so they are formed by the stage iteration rather than here.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    logical, intent(in) :: prevs
        !! Defines the status of the previous step.  The value is true
        !! if the previous step was successful; else, false if the
        !! previous step failed.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives
        !! at x.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Process
    return
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_post_step(this, sys, dense, x, xn, y, yn, f, fn, k, args)
    !! Retains the first-stage derivative for cubic Hermite interpolation.
    !!
    !! For a mass-matrix problem, the driver's f is a raw RHS rather than a
    !! derivative. The first stage supplies the consistent left-end slope;
    !! the stiffly accurate final stage already supplies the right-end slope.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    logical, intent(in) :: dense
        !! Determines if dense output is requested (true); else, false.
    real(real64), intent(in) :: x
        !! The previous value of the independent variable.
    real(real64), intent(in) :: xn
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the solution at x.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the derivatives at x.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.
    real(real64), intent(inout), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    this%m_startDerivative = k(:,1)
end subroutine

! ------------------------------------------------------------------------------
pure function kc_get_is_fsal(this) result(rst)
    !! Gets a logical parameter stating if this is a first-same-as-last
    !! (FSAL) integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    logical :: rst
        !! True for a FSAL integrator; else, false.
    rst = .false.
end function

! ------------------------------------------------------------------------------
pure function kc_get_stage_count(this) result(rst)
    !! Gets the stage count for this integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    integer :: rst
        !! The stage count.
    select type (this)
    type is (kennedy_carpenter_4)
        rst = 6
    type is (kennedy_carpenter_5)
        rst = 8
    class default
        rst = 0
    end select
end function

! ------------------------------------------------------------------------------
pure function kc4_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(kennedy_carpenter_4), intent(in) :: this
        !! The kennedy_carpenter_4 object.
    integer :: rst
        !! The order.
    rst = 4
end function

! ------------------------------------------------------------------------------
pure function kc5_get_order(this) result(rst)
    !! Gets the order of the integrator.
    class(kennedy_carpenter_5), intent(in) :: this
        !! The kennedy_carpenter_5 object.
    integer :: rst
        !! The order.
    rst = 5
end function

! ------------------------------------------------------------------------------
subroutine kc_attempt_step(this, sys, h, x, y, f, yn, fn, yerr, k, args)
    !! Attempts an integration step for this integrator.
    !!
    !! The state at each stage is solved from
    !! \f[
    !! M(Y_i-Y_i^*)-h a_{ii}f(x+c_i h,Y_i)=0,
    !! \f]
    !! using Newton iteration.  The high- and embedded-order solutions are
    !! \f$y_{n+1}=y_n+h\sum b_i k_i\f$ and
    !! \f$\widehat y_{n+1}=y_n+h\sum d_i k_i\f$, respectively.
    class(kennedy_carpenter), intent(inout) :: this
        !! The kennedy_carpenter object.
    class(ode_container), intent(inout) :: sys
        !! The ode_container object containing the ODE's to integrate.
    real(real64), intent(in) :: h
        !! The current step size.
    real(real64), intent(in) :: x
        !! The current value of the independent variable.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the current solution at x.
    real(real64), intent(in), dimension(:) :: f
        !! An N-element array containing the values of the derivatives
        !! at x.
    real(real64), intent(out), dimension(:) :: yn
        !! An N-element array where this routine will write the next
        !! solution estimate at x + h.
    real(real64), intent(out), dimension(:) :: fn
        !! An N-element array where this routine will write the next
        !! derivative estimate at x + h.
    real(real64), intent(out), dimension(:) :: yerr
        !! An N-element array where this routine will write an estimate
        !! of the error in each equation.
    real(real64), intent(out), dimension(:,:) :: k
        !! An N-by-NSTAGES matrix containing the derivatives at each stage.
    class(*), intent(inout), optional :: args
        !! An optional argument that can be used to pass information
        !! in and out of the differential equation subroutine.

    ! Local Variables
    integer :: i, j, n, stages, njac
    logical :: usemass, massDependent
    real(real64) :: a(8,8), b(8), d(8), c(8)
    real(real64), allocatable, dimension(:) :: base, state, embedded

    ! Initialization
    n = size(y)
    usemass = sys%has_mass_matrix()
    massDependent = sys%get_is_mass_matrix_dependent()
    call kc_table(this, a, b, d, c, stages)
    allocate( &
        base(n), &
        state(n), &
        embedded(n) &
    )

    ! Process
    ! The first stage of an ESDIRK method is explicit
    k = 0.0d0
    if (usemass) then
        if (allocated(this%m_mass)) then
            if (size(this%m_mass, 1) /= n) then
                deallocate(this%m_mass)
                this%m_massComputed = .false.
            end if
        end if
        if (.not.allocated(this%m_mass)) allocate(this%m_mass(n, n))
        if (.not.this%m_massComputed .or. massDependent) then
            call sys%compute_mass_matrix(x, y, this%m_mass, args)
            this%m_massComputed = .true.
        end if
        call kc_consistent_derivative(sys, x, y, f, this%m_mass, k(:,1), args)
        call this%log_jacobian_evaluation()
        call this%log_function_evaluation()
        call this%log_matrix_factorization()
        if (sys%use_finite_difference_jacobian()) then
            call this%log_function_evaluation(n + 1)
        end if
    else
        k(:,1) = f
    end if

    ! Each remaining stage is diagonally implicit
    do i = 2, stages
        ! Accumulate the contributions from the previously computed stages
        base = y
        do j = 1, i - 1
            base = base + h * a(i,j) * k(:,j)
        end do

        if (usemass) then
            call kc_mass_stage(sys, x + c(i) * h, base, h * a(i,i), &
                this%m_mass, massDependent, state, k(:,i), njac, args)
            call this%log_function_evaluation(njac + 1)
        else
            call kc_ode_stage(sys, x + c(i) * h, base, h * a(i,i), &
                state, k(:,i), njac, args)
            call this%log_function_evaluation(njac + 2)
        end if
        call this%log_jacobian_evaluation(njac)
        call this%log_matrix_factorization(njac)
        if (sys%use_finite_difference_jacobian()) then
            call this%log_function_evaluation(njac * (n + 1))
        end if
    end do

    ! Form the solution estimate and the embedded solution estimate
    yn = y
    embedded = y
    do i = 1, stages
        yn = yn + h * b(i) * k(:,i)
        embedded = embedded + h * d(i) * k(:,i)
    end do
    yerr = yn - embedded

    ! Both tableaus are stiffly accurate.  The final stage derivative is
    ! therefore the derivative at the accepted state, even when M is singular.
    fn = k(:, stages)
end subroutine

! ------------------------------------------------------------------------------
subroutine kc_consistent_derivative(sys, x, y, f, mass, derivative, args)
    !! Computes a consistent initial derivative for an index-1 DAE.
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, y(:), f(:), mass(:,:)
    real(real64), intent(out) :: derivative(:)
    class(*), intent(inout), optional :: args
    real(real64), allocatable :: singular_values(:), left_vectors(:,:), &
        jac(:,:), fplus(:), augmented(:,:), rhs(:)
    real(real64) :: fdstep, rank_tolerance, scale
    integer :: n, rank, nullity, i

    n = size(y)
    call svd(mass, singular_values, left_vectors)
    scale = max(1.0d0, singular_values(1))
    rank_tolerance = 100.0d0 * epsilon(1.0d0) * scale
    rank = count(singular_values > rank_tolerance)

    if (rank == n) then
        derivative = solve_kc_system(mass, f)
        return
    end if

    nullity = n - rank
    allocate(jac(n,n), fplus(n), augmented(n,n), rhs(n))
    call sys%compute_jacobian(x, y, jac, args)
    fdstep = sys%get_finite_difference_step()
    call sys%compute_rhs(x + fdstep, y, fplus, args)
    fplus = (fplus - f) / fdstep

    augmented = 0.0d0
    rhs = 0.0d0
    augmented(1:rank,:) = matmul(transpose(left_vectors(:,1:rank)), mass)
    rhs(1:rank) = matmul(transpose(left_vectors(:,1:rank)), f)
    do i = 1, nullity
        augmented(rank+i,:) = matmul(left_vectors(:,rank+i), jac)
        rhs(rank+i) = -dot_product(left_vectors(:,rank+i), fplus)
    end do

    derivative = solve_kc_system(augmented, rhs)
end subroutine

subroutine kc_ode_stage(sys, x, base, diagonal_step, state, derivative, njac, args)
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, diagonal_step
    real(real64), intent(in), dimension(:) :: base
    real(real64), intent(out), dimension(:) :: state, derivative
    integer, intent(out) :: njac
    class(*), intent(inout), optional :: args
    real(real64), dimension(size(base)) :: rhs, residual, delta
    real(real64), dimension(size(base),size(base)) :: jac, system
    integer :: equation, iteration

    call sys%compute_rhs(x, base, rhs, args)
    state = base + diagonal_step * rhs
    njac = 0
    do iteration = 1, 12
        call sys%compute_rhs(x, state, rhs, args)
        residual = state - base - diagonal_step * rhs
        if (norm2(residual) <= 1.0d-12 * max(1.0d0, norm2(state))) exit
        call sys%compute_jacobian(x, state, jac, args)
        njac = njac + 1
        system = -diagonal_step * jac
        do equation = 1, size(base)
            system(equation,equation) = system(equation,equation) + 1.0d0
        end do
        delta = solve_kc_system(system, -residual)
        state = state + delta
    end do
    if (iteration > 12) error stop DIFFEQ_CONVERGENCE_ERROR
    derivative = rhs
end subroutine

subroutine kc_mass_stage(sys, x, base, diagonal_step, mass, massDependent, &
    state, derivative, njac, args)
    !! Solves one Kennedy--Carpenter mass-matrix stage without forming M^{-1}f.
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, base(:), diagonal_step
    real(real64), intent(inout), dimension(:,:) :: mass
    logical, intent(in) :: massDependent
    real(real64), intent(out) :: state(:), derivative(:)
    integer, intent(out) :: njac
    class(*), intent(inout), optional :: args
    real(real64) :: f(size(base)), residual(size(base))
    real(real64) :: jac(size(base),size(base))
    real(real64) :: mass_perturbed(size(base),size(base))
    real(real64) :: system(size(base),size(base)), delta(size(base))
    real(real64) :: fdstep
    integer :: i, iteration

    state = base
    njac = 0
    do iteration = 1, 12
        derivative = (state - base) / diagonal_step
        call sys%compute_rhs(x, state, f, args)
        if (massDependent) call sys%compute_mass_matrix(x, state, mass, args)
        residual = matmul(mass, state - base) - diagonal_step * f
        if (norm2(residual) <= 1.0d-14 * max(1.0d0, norm2(state))) exit
        call sys%compute_jacobian(x, state, jac, args)
        njac = njac + 1
        system = mass - diagonal_step * jac
        if (massDependent) then
            fdstep = sys%get_finite_difference_step()
            do i = 1, size(base)
                state(i) = state(i) + fdstep
                call sys%compute_mass_matrix(x, state, mass_perturbed, args)
                state(i) = state(i) - fdstep
                system(:,i) = system(:,i) + &
                    matmul((mass_perturbed - mass) / fdstep, state - base)
            end do
        end if
        delta = solve_kc_system(system, -residual)
        state = state + delta
    end do
    if (iteration > 12) error stop DIFFEQ_CONVERGENCE_ERROR
    derivative = (state - base) / diagonal_step
end subroutine

subroutine kc_table(this, a, b, d, c, stages)
    !! Populates the Butcher tableau for the requested integrator.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(out) :: a(8,8)
        !! The matrix of stage coefficients.  Only the leading
        !! stages-by-stages block is populated.
    real(real64), intent(out) :: b(8)
        !! The array of weights defining the higher-order solution.
    real(real64), intent(out) :: d(8)
        !! The array of weights defining the embedded solution.
    real(real64), intent(out) :: c(8)
        !! The array of nodes at which each stage is evaluated.
    integer, intent(out) :: stages
        !! The number of stages used by the integrator.

    ! Initialization
    a = 0.0d0
    b = 0.0d0
    d = 0.0d0
    c = 0.0d0

    ! Process
    select type (this)
    type is (kennedy_carpenter_4)
        ! ARK4(3)6L[2]SA
        stages = 6

        a(2,1) = 1.0d0 / 4.0d0
        a(2,2) = 1.0d0 / 4.0d0

        a(3,1) = 8611.0d0 / 62500.0d0
        a(3,2) = -1743.0d0 / 31250.0d0
        a(3,3) = 1.0d0 / 4.0d0

        a(4,1) = 5012029.0d0 / 34652500.0d0
        a(4,2) = -654441.0d0 / 2922500.0d0
        a(4,3) = 174375.0d0 / 388108.0d0
        a(4,4) = 1.0d0 / 4.0d0

        a(5,1) = 15267082809.0d0 / 155376265600.0d0
        a(5,2) = -71443401.0d0 / 120774400.0d0
        a(5,3) = 730878875.0d0 / 902184768.0d0
        a(5,4) = 2285395.0d0 / 8070912.0d0
        a(5,5) = 1.0d0 / 4.0d0

        a(6,1) = 82889.0d0 / 524892.0d0
        a(6,3) = 15625.0d0 / 83664.0d0
        a(6,4) = 69875.0d0 / 102672.0d0
        a(6,5) = -2260.0d0 / 8211.0d0
        a(6,6) = 1.0d0 / 4.0d0

        ! The method is stiffly accurate, so the weights match the last stage
        b(1) = a(6,1)
        b(3) = a(6,3)
        b(4) = a(6,4)
        b(5) = a(6,5)
        b(6) = a(6,6)

        d(1) = 4586570599.0d0 / 29645900160.0d0
        d(3) = 178811875.0d0 / 945068544.0d0
        d(4) = 814220225.0d0 / 1159782912.0d0
        d(5) = -3700637.0d0 / 11593932.0d0
        d(6) = 61727.0d0 / 225920.0d0

        c(2) = 0.5d0
        c(3) = 83.0d0 / 250.0d0
        c(4) = 31.0d0 / 50.0d0
        c(5) = 17.0d0 / 20.0d0
        c(6) = 1.0d0
    type is (kennedy_carpenter_5)
        ! ARK5(4)8L[2]SA
        stages = 8

        a(2,1) = 41.0d0 / 200.0d0
        a(2,2) = 41.0d0 / 200.0d0

        a(3,1) = 41.0d0 / 400.0d0
        a(3,2) = -567603406766.0d0 / 11931857230679.0d0
        a(3,3) = 41.0d0 / 200.0d0

        a(4,1) = 683785636431.0d0 / 9252920307686.0d0
        a(4,3) = -110385047103.0d0 / 1367015193373.0d0
        a(4,4) = 41.0d0 / 200.0d0

        a(5,1) = 3016520224154.0d0 / 10081342136671.0d0
        a(5,3) = 30586259806659.0d0 / 12414158314087.0d0
        a(5,4) = -22760509404356.0d0 / 11113319521817.0d0
        a(5,5) = 41.0d0 / 200.0d0

        a(6,1) = 218866479029.0d0 / 1489978393911.0d0
        a(6,3) = 638256894668.0d0 / 5436446318841.0d0
        a(6,4) = -1179710474555.0d0 / 5321154724896.0d0
        a(6,5) = -60928119172.0d0 / 8023461067671.0d0
        a(6,6) = 41.0d0 / 200.0d0

        a(7,1) = 1020004230633.0d0 / 5715676835656.0d0
        a(7,3) = 25762820946817.0d0 / 25263940353407.0d0
        a(7,4) = -2161375909145.0d0 / 9755907335909.0d0
        a(7,5) = -211217309593.0d0 / 5846859502534.0d0
        a(7,6) = -4269925059573.0d0 / 7827059040749.0d0
        a(7,7) = 41.0d0 / 200.0d0

        a(8,1) = -872700587467.0d0 / 9133579230613.0d0
        a(8,4) = 22348218063261.0d0 / 9555858737531.0d0
        a(8,5) = -1143369518992.0d0 / 8141816002931.0d0
        a(8,6) = -39379526789629.0d0 / 19018526304540.0d0
        a(8,7) = 32727382324388.0d0 / 42900044865799.0d0
        a(8,8) = 41.0d0 / 200.0d0

        ! The method is stiffly accurate, so the weights match the last stage
        b(1) = a(8,1)
        b(4) = a(8,4)
        b(5) = a(8,5)
        b(6) = a(8,6)
        b(7) = a(8,7)
        b(8) = a(8,8)

        d(1) = -975461918565.0d0 / 9796059967033.0d0
        d(4) = 78070527104295.0d0 / 32432590147079.0d0
        d(5) = -548382580838.0d0 / 3424219808633.0d0
        d(6) = -33438840321285.0d0 / 15594753105479.0d0
        d(7) = 3629800801594.0d0 / 4656183773603.0d0
        d(8) = 4035322873751.0d0 / 18575991585200.0d0

        c(2) = 41.0d0 / 100.0d0
        c(3) = 2935347310677.0d0 / 11292855782101.0d0
        c(4) = 1426016391358.0d0 / 7196633302097.0d0
        c(5) = 0.92d0
        c(6) = 0.24d0
        c(7) = 0.6d0
        c(8) = 1.0d0
    end select
end subroutine

! ------------------------------------------------------------------------------
function solve_kc_system(a, b) result(x)
    !! Solves the linear system \f$ A x = b \f$ by means of a QR factorization
    !! with column pivoting.
    real(real64), intent(in), dimension(:,:) :: a
        !! The N-by-N system matrix.
    real(real64), intent(in), dimension(:) :: b
        !! An N-element array containing the right-hand side.
    real(real64), allocatable, dimension(:) :: x
        !! An N-element array containing the solution.

    ! Process
    x = solve_least_squares_full(a, b)
end function

! ------------------------------------------------------------------------------
subroutine kc_interpolate(this, x, xn, yn, fn, xn1, yn1, fn1, y)
    !! Performs the interpolation.
    !!
    !! A cubic Hermite polynomial is constructed from the solution and
    !! derivative values at each end of the step.
    class(kennedy_carpenter), intent(in) :: this
        !! The kennedy_carpenter object.
    real(real64), intent(in) :: x
        !! The value of the independent variable at which to compute
        !! the interpolation.
    real(real64), intent(in) :: xn
        !! The previous value of the independent variable at which the
        !! solution is computed.
    real(real64), intent(in), dimension(:) :: yn
        !! An N-element array containing the solution at xn.
    real(real64), intent(in), dimension(:) :: fn
        !! An N-element array containing the derivatives at xn.
    real(real64), intent(in) :: xn1
        !! The value of the independent variable at xn + h.
    real(real64), intent(in), dimension(:) :: yn1
        !! An N-element array containing the solution at xn + h.
    real(real64), intent(in), dimension(:) :: fn1
        !! An N-element array containing the derivatives at xn + h.
    real(real64), intent(out), dimension(:) :: y
        !! An N-element array where this routine will write the
        !! solution values interpolated at x.

    ! Local Variables
    real(real64) :: h, s

    ! Initialization
    h = xn1 - xn
    s = (x - xn) / h

    ! Process
    y = (2.0d0 * s**3 - 3.0d0 * s**2 + 1.0d0) * yn + &
        (-2.0d0 * s**3 + 3.0d0 * s**2) * yn1 + &
        (s**3 - s**2) * h * fn1
    if (allocated(this%m_startDerivative)) then
        y = y + (s**3 - 2.0d0 * s**2 + s) * h * this%m_startDerivative
    else
        y = y + (s**3 - 2.0d0 * s**2 + s) * h * fn
    end if
end subroutine

! ------------------------------------------------------------------------------
end module