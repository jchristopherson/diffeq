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

module diffeq_index_reduction
    !! Stabilized index-1 Gear-Gupta-Leimkuhler-inspired (GGL) reduction
    !! for mechanical systems with autonomous holonomic constraints.
    !!
    !! For autonomous holonomic constraints \(g(q)=0\), this variant uses
    !! coordinates \(q\), velocities \(v\), force multipliers \(\lambda\),
    !! and kinematic multipliers \(\mu\) in the formulation
    !! \[
    !! \dot q = v + G(q)^T\mu, \qquad
    !! M(q)\dot v = F(x,q,v) - G(q)^T\lambda,
    !! \]
    !! \[
    !! 0 = G\dot q + \alpha g, \qquad
    !! 0 = DG(q)[\dot q]v + GM^{-1}(F-G^T\lambda) + \beta Gv,
    !! \]
    !! where \(G=\partial g/\partial q\). Thus \(\dot g=-\alpha g\)
    !! and \(d(Gv)/dx=-\beta Gv\). Consistent constraints remain zero
    !! in the exact solution; positive gains damp violations. Numerical
    !! integration and interpolation do not project onto the constraints.
    !!
    !! Unlike classical index-2 GGL, the original position and velocity
    !! constraints are replaced by their stabilized evolution equations.
    !! The algebraic multipliers are locally determined when \(GG^T\)
    !! and \(GM^{-1}G^T\) are nonsingular. Independent constraints and a
    !! symmetric positive-definite mechanical mass matrix suffice. Explicitly
    !! moving constraints are not supported: the constraint callbacks must
    !! not depend on the independent variable.
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use diffeq_base
    use diffeq_errors
    use linalg, only : identity, lu_factor, solve_lu
    implicit none
    private
    public :: ggl_container
    public :: dae_algebraic_constraints
    public :: dae_constraint_curvature

    abstract interface
        subroutine dae_algebraic_constraints(x, y, g, args)
            !! Evaluates algebraic constraint residuals \(g(x,y)\).
            !!
            !! For \(N\) input state components and \(M\) constraints,
            !! the callback represents a mapping
            !! \[
            !! g : \mathbb{R}\times\mathbb{R}^{N}
            !! \longrightarrow \mathbb{R}^{M}.
            !! \]
            !! A consistent state satisfies \(g(x,y)=0\). Return the
            !! residuals themselves, not their derivatives, a Boolean result,
            !! or a projected state. A nonzero residual is a valid evaluation
            !! at an inconsistent state, rather than a callback failure.
            !!
            !! The caller supplies the input and output array sizes. The
            !! implementation must interpret the same state ordering as its
            !! model and write every output element. This interface performs
            !! no dimension, rank, or initial-consistency checks. Residual
            !! ordering and physical units are defined by the model; any
            !! subsequent consistency tolerance must respect that scaling.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: x
                !! The independent-variable value at which to evaluate the
                !! constraints, normally time. Autonomous constraints may
                !! ignore this argument.
            real(real64), intent(in), dimension(:) :: y
                !! An N-element input state in the ordering defined by the
                !! model. No particular coordinate, velocity, or multiplier
                !! partition is imposed by this interface. It must not be
                !! modified by the callback.
            real(real64), intent(out), dimension(:) :: g
                !! An M-element output array receiving the constraint
                !! residuals in the model's declared constraint order.
                !! Every element must be assigned on every call; incoming
                !! contents are undefined. Exact constraint satisfaction
                !! corresponds to a zero vector.
            class(*), intent(inout), optional :: args
                !! Optional user-owned model parameters or evaluation state.
                !! Check presence and dynamic type before accessing it.
                !! The callback may update this object, but repeated calls
                !! with the same inputs should give the same residuals unless
                !! the user intentionally changes model parameters. This
                !! argument is not reserved for index-reduction workspace.
        end subroutine

        subroutine dae_constraint_curvature(x, q, direction, v, curvature, args)
            !! Evaluates \(DG(q)[\mathrm{direction}]v\), without including
            !! \(G\dot v\). Each component is
            !! \(\mathrm{direction}^T\nabla^2 g_i(q)v\).
            !! Assign every output element. For linear constraints return
            !! zeros. The callback must support off-manifold coordinates and
            !! arbitrary directions used during nonlinear/Jacobian evaluation.
            import real64
            real(real64), intent(in) :: x
                !! Independent-variable value, provided for callback uniformity;
                !! autonomous constraints must not depend on it.
            real(real64), intent(in), dimension(:) :: q
                !! N coordinates, in the same order as the force model.
            real(real64), intent(in), dimension(:) :: direction
                !! N components of the coordinate derivative \(\dot q\),
                !! including the kinematic multiplier correction.
            real(real64), intent(in), dimension(:) :: v
                !! N physical velocity components; these need not equal direction.
            real(real64), intent(out), dimension(:) :: curvature
                !! M curvature terms in the order of the constraint callback.
            class(*), intent(inout), optional :: args
                !! Optional user-owned parameters, forwarded unchanged.
        end subroutine
    end interface

    type, extends(ode_container) :: ggl_container
        !! Mechanical model adapter implementing stabilized index-1 GGL.
        !!
        !! Use this type instead of [[ode_container]] when reducing an index-3
        !! holonomically constrained mechanical model. Call `initialize` and
        !! assign all required callbacks before using an index-1 DAE solver.
        !! For N coordinates and M constraints, the state and RHS layout is
        !! \[
        !! z=[q_N,v_N,\lambda_M,\mu_M]^T, \quad
        !! \mathcal M=\operatorname{diag}(I_N,M_N,0_M,0_M),
        !! \]
        !! \[
        !! \mathcal f=\begin{bmatrix}
        !! v+G^T\mu \\
        !! F-G^T\lambda \\
        !! G(v+G^T\mu)+\alpha g \\
        !! DG(q)[\dot q]v+GM^{-1}(F-G^T\lambda)+\beta Gv
        !! \end{bmatrix}.
        !! \]
        !! All constraint vectors use the same M-component ordering.
        !!
        !! The inherited `fcn` has a specialized contract here: its input is
        !! [q,v] of length 2N and its output is the applied force F of length N,
        !! excluding multiplier forces. The optional inherited `mass_matrix`
        !! receives q of length N and writes the physical N-by-N mass matrix;
        !! if absent, the physical matrix is identity. Set the inherited mass
        !! dependency flag true whenever this matrix varies with x or q.
        !!
        !! An optional inherited `jacobian` must differentiate the complete
        !! augmented RHS with respect to z, not just the applied force. Without
        !! it, the inherited finite-difference implementation differentiates
        !! this container's `compute_rhs`. User args are never used as internal
        !! workspace. Callbacks must accept perturbed, off-manifold states.
        integer, private :: m_coordinates = 0
        integer, private :: m_constraints = 0
        real(real64), private :: m_alpha = 1.0d0
        real(real64), private :: m_beta = 1.0d0
        procedure(dae_algebraic_constraints), pointer, public, nopass :: constraints => null()
            !! Required callback receiving N coordinates q and writing M
            !! autonomous position residuals g(q), including off-manifold values.
        procedure(ode_jacobian), pointer, public, nopass :: constraint_jacobian => null()
            !! Required callback receiving q and writing the M-by-N matrix
            !! \(G=\partial g/\partial q\). No finite-difference fallback is used.
        procedure(dae_constraint_curvature), pointer, public, nopass :: constraint_curvature => null()
            !! Required callback for \(DG(q)[\dot q]v\).
    contains
        procedure :: initialize => ggl_initialize
        procedure :: set_stabilization_parameters => ggl_set_stabilization
        procedure :: get_state_count => ggl_state_count
        procedure :: has_mass_matrix => ggl_has_mass
        procedure :: get_is_ode_defined => ggl_is_defined
        procedure :: compute_rhs => ggl_rhs
        procedure :: compute_mass_matrix => ggl_mass
        procedure :: validate_initial_conditions => ggl_validate
        procedure :: initialize_multipliers => ggl_initialize_multipliers
        procedure, private :: check_state => ggl_check_state
        procedure, private :: physical_mass => ggl_physical_mass
    end type

contains
subroutine ggl_initialize(this, coordinate_count, constraint_count, alpha, beta)
    !! Defines the mechanical dimensions and resets stabilization gains to one
    !! unless supplied. Callback associations and inherited model options are
    !! retained. Invalid counts or gains cause DIFFEQ_INVALID_INPUT_ERROR.
    class(ggl_container), intent(inout) :: this
        !! Container to configure.
    integer, intent(in) :: coordinate_count
        !! Positive number N of coordinates and physical velocities.
    integer, intent(in) :: constraint_count
        !! Number M of independent constraints, with 1 <= M <= N.
    real(real64), intent(in), optional :: alpha, beta
        !! Nonnegative finite position and velocity decay rates, respectively,
        !! in inverse independent-variable units. Zero disables the respective
        !! damping term without disabling the differentiated constraint.
    real(real64) :: position_rate, velocity_rate

    if (coordinate_count < 1 .or. constraint_count < 1 .or. &
        constraint_count > coordinate_count) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_coordinates = coordinate_count
    this%m_constraints = constraint_count
    position_rate = 1.0d0
    velocity_rate = 1.0d0
    if (present(alpha)) position_rate = alpha
    if (present(beta)) velocity_rate = beta
    call this%set_stabilization_parameters(position_rate, velocity_rate)
end subroutine

subroutine ggl_set_stabilization(this, alpha, beta)
    !! Sets the two independent constraint decay rates. Negative or nonfinite
    !! values cause DIFFEQ_INVALID_INPUT_ERROR. No state projection is performed.
    class(ggl_container), intent(inout) :: this
        !! Configured or unconfigured container whose gains are to be changed.
    real(real64), intent(in) :: alpha, beta
        !! Position and velocity decay rates, respectively; both have units
        !! of inverse time when x is time. Zero preserves rather than damps
        !! the corresponding residual in the exact solution.

    if (.not.ieee_is_finite(alpha) .or. .not.ieee_is_finite(beta)) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    if (alpha < 0.0d0 .or. beta < 0.0d0) error stop DIFFEQ_INVALID_INPUT_ERROR
    this%m_alpha = alpha
    this%m_beta = beta
end subroutine

function ggl_state_count(this) result(rst)
    !! Returns 2N+2M, or zero before dimensions have been initialized.
    class(ggl_container), intent(in) :: this
        !! Container whose augmented state size is queried.
    integer :: rst
        !! Required length of the state, initial conditions, and RHS arrays.

    rst = 2 * (this%m_coordinates + this%m_constraints)
end function

function ggl_has_mass(this) result(rst)
    !! Always returns true: the augmented matrix is singular even when the
    !! physical matrix defaults to identity. Explicit ODE paths are unsuitable.
    class(ggl_container), intent(in) :: this
        !! Container being queried; initialization is not required.
    logical :: rst
        !! True, indicating that a mass-matrix DAE path is required.

    rst = .true.
end function

function ggl_is_defined(this) result(rst)
    !! Checks configuration and all required callback associations. This does
    !! not evaluate callbacks or certify constraint rank or mass regularity.
    class(ggl_container), intent(in) :: this
        !! Container whose model definition is inspected.
    logical :: rst
        !! True only for initialized dimensions and associated force, position,
        !! constraint-Jacobian, and curvature callbacks. Mass is optional.

    rst = this%m_coordinates > 0 .and. this%m_constraints > 0 .and. &
        associated(this%fcn) .and. associated(this%constraints) .and. &
        associated(this%constraint_jacobian) .and. associated(this%constraint_curvature)
end function

subroutine ggl_check_state(this, y)
    !! Checks dimensions and finite state values without invoking callbacks.
    !! Uninitialized dimensions/nonfinite values cause INVALID_INPUT; a wrong
    !! array length causes DIFFEQ_ARRAY_SIZE_ERROR.
    class(ggl_container), intent(in) :: this
        !! Container defining the required state dimensions.
    real(real64), intent(in), dimension(:) :: y
        !! Augmented state [q,v,lambda,mu] of length 2N+2M.

    if (this%m_coordinates < 1 .or. this%m_constraints < 1) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    if (size(y) /= this%get_state_count()) error stop DIFFEQ_ARRAY_SIZE_ERROR
    if (.not.all(ieee_is_finite(y))) error stop DIFFEQ_INVALID_INPUT_ERROR
end subroutine

subroutine ggl_physical_mass(this, x, q, mass, args)
    !! Evaluates the physical mass callback, or supplies identity. This private
    !! helper assumes the caller has checked array dimensions. Singular masses
    !! are diagnosed by the linear solve, not by this evaluation.
    class(ggl_container), intent(inout) :: this
        !! Container holding the optional mechanical mass callback.
    real(real64), intent(in) :: x
        !! Independent-variable value for the mechanical mass callback.
    real(real64), intent(in), dimension(:) :: q
        !! N coordinates, without velocities or multipliers.
    real(real64), intent(out), dimension(:,:) :: mass
        !! Physical N-by-N matrix, not the augmented singular matrix.
    class(*), intent(inout), optional :: args
        !! Optional user-owned parameters forwarded to the callback.

    if (associated(this%mass_matrix)) then
        call this%mass_matrix(x, q, mass, args)
    else
        mass = identity(size(q))
    end if
end subroutine

subroutine ggl_rhs(this, x, y, dydx, args)
    !! Assembles the four RHS blocks documented on [[ggl_container]]. The
    !! output is a mass-matrix RHS, not the full derivative z'. Only the
    !! physical mass is resolved to compute the acceleration residual.
    !! Missing callbacks cause NULL_POINTER; invalid state/output lengths
    !! cause INVALID_INPUT/ARRAY_SIZE; singular physical mass causes SINGULAR_MATRIX.
    class(ggl_container), intent(inout) :: this
        !! Initialized container with all required callbacks assigned.
    real(real64), intent(in) :: x
        !! Current independent-variable value.
    real(real64), intent(in), dimension(:) :: y
        !! State [q,v,lambda,mu] of length 2N+2M; off-manifold values are allowed.
    real(real64), intent(out), dimension(:) :: dydx
        !! 2N+2M RHS entries, ordered as kinematics, momentum, stabilized
        !! position residual, and stabilized velocity residual.
    class(*), intent(inout), optional :: args
        !! Optional user parameters forwarded unchanged to every callback.
    integer :: coordinates, constraints, physical
    real(real64), allocatable, dimension(:) :: position, curvature, qdot
    real(real64), allocatable, dimension(:,:) :: gradient, mass, acceleration

    call this%check_state(y)
    if (.not.this%get_is_ode_defined()) error stop DIFFEQ_NULL_POINTER_ERROR
    if (size(dydx) /= size(y)) error stop DIFFEQ_ARRAY_SIZE_ERROR
    coordinates = this%m_coordinates
    constraints = this%m_constraints
    physical = 2 * coordinates
    allocate(position(constraints), curvature(constraints), qdot(coordinates), &
        gradient(constraints,coordinates), mass(coordinates,coordinates), &
        acceleration(coordinates,1))
    call this%constraints(x, y(1:coordinates), position, args)
    call this%constraint_jacobian(x, y(1:coordinates), gradient, args)
    qdot = y(coordinates+1:physical) + &
        matmul(transpose(gradient), y(physical+constraints+1:))
    dydx(1:coordinates) = qdot
    call this%fcn(x, y(1:physical), dydx(coordinates+1:physical), args)
    dydx(coordinates+1:physical) = dydx(coordinates+1:physical) - &
        matmul(transpose(gradient), y(physical+1:physical+constraints))
    call this%physical_mass(x, y(1:coordinates), mass, args)
    acceleration(:,1) = dydx(coordinates+1:physical)
    call ggl_linear_solve(mass, acceleration)
    call this%constraint_curvature(x, y(1:coordinates), qdot, &
        y(coordinates+1:physical), curvature, args)
    dydx(physical+1:physical+constraints) = matmul(gradient, qdot) + this%m_alpha * position
    dydx(physical+constraints+1:) = curvature + matmul(gradient, acceleration(:,1)) + &
        this%m_beta * matmul(gradient, y(coordinates+1:physical))
end subroutine

subroutine ggl_mass(this, x, y, m, args)
    !! Assembles diag(I,M,0,0). Both multiplier blocks and their cross blocks
    !! are zero. No constraints or force callbacks are needed for this operation.
    !! Uninitialized/wrong state sizes or matrix sizes are diagnosed explicitly.
    class(ggl_container), intent(inout) :: this
        !! Container defining dimensions and the optional physical mass.
    real(real64), intent(in) :: x
        !! Current independent-variable value.
    real(real64), intent(in), dimension(:) :: y
        !! Augmented state [q,v,lambda,mu] of length 2N+2M.
    real(real64), intent(out), dimension(:,:) :: m
        !! Square augmented matrix with side length 2N+2M. A wrong shape
        !! causes DIFFEQ_MATRIX_SIZE_ERROR.
    class(*), intent(inout), optional :: args
        !! Optional user parameters forwarded to the physical mass callback.
    integer :: coordinates

    call this%check_state(y)
    if (size(m,1) /= size(y) .or. size(m,2) /= size(y)) &
        error stop DIFFEQ_MATRIX_SIZE_ERROR
    coordinates = this%m_coordinates
    m = 0.0d0
    m(1:coordinates,1:coordinates) = identity(coordinates)
    call this%physical_mass(x, y(1:coordinates), &
        m(coordinates+1:2*coordinates,coordinates+1:2*coordinates), args)
end subroutine

function ggl_validate(this, x, y, tolerance, args) result(rst)
    !! Checks max(abs(g)) and max(abs(Gv)) against one absolute tolerance.
    !! This tests physical position/velocity consistency, not the multiplier
    !! equations, and never changes the state. Choose units/scaling so one
    !! tolerance is meaningful for both residuals. Nonfinite residuals fail.
    class(ggl_container), intent(inout) :: this
        !! Initialized container; only constraints and their Jacobian are needed.
    real(real64), intent(in) :: x
        !! Independent-variable value passed to the callbacks.
    real(real64), intent(in), dimension(:) :: y
        !! Augmented initial state of length 2N+2M.
    real(real64), intent(in) :: tolerance
        !! Finite nonnegative absolute residual threshold; invalid values
        !! cause DIFFEQ_INVALID_INPUT_ERROR.
    class(*), intent(inout), optional :: args
        !! Optional user parameters forwarded without using internal workspace.
    logical :: rst
        !! True when both physical residuals are finite and within tolerance.
    real(real64), allocatable, dimension(:) :: position, velocity
    real(real64), allocatable, dimension(:,:) :: gradient

    call this%check_state(y)
    if (.not.ieee_is_finite(tolerance)) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (tolerance < 0.0d0) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.associated(this%constraints) .or. .not.associated(this%constraint_jacobian)) &
        error stop DIFFEQ_NULL_POINTER_ERROR
    allocate(position(this%m_constraints), velocity(this%m_constraints), &
        gradient(this%m_constraints,this%m_coordinates))
    call this%constraints(x, y(1:this%m_coordinates), position, args)
    call this%constraint_jacobian(x, y(1:this%m_coordinates), gradient, args)
    velocity = matmul(gradient, y(this%m_coordinates+1:2*this%m_coordinates))
    rst = all(ieee_is_finite(position)) .and. all(ieee_is_finite(velocity))
    rst = rst .and. maxval(abs(position)) <= tolerance .and. maxval(abs(velocity)) <= tolerance
end function

subroutine ggl_initialize_multipliers(this, x, y, args)
    !! Computes mu and lambda to satisfy the two stabilized algebraic rows
    !! for the supplied q and v. Coordinates and physical velocities are not
    !! projected or changed, even if they violate the original constraints.
    !! Singular GG^T, physical mass, or GM^{-1}G^T causes SINGULAR_MATRIX.
    class(ggl_container), intent(inout) :: this
        !! Initialized container with all required callbacks assigned.
    real(real64), intent(in) :: x
        !! Initial independent-variable value.
    real(real64), intent(inout), dimension(:) :: y
        !! Full [q,v,lambda,mu] state. Input multipliers are ignored and replaced;
        !! output multipliers are algebraically consistent to solve accuracy.
    class(*), intent(inout), optional :: args
        !! Optional user-owned parameters forwarded to model callbacks.
    integer :: coordinates, constraints, physical
    real(real64), allocatable, dimension(:) :: position, rhs
    real(real64), allocatable, dimension(:,:) :: gradient, mass, inverse_gradient, multiplier

    call this%check_state(y)
    if (.not.this%get_is_ode_defined()) error stop DIFFEQ_NULL_POINTER_ERROR
    coordinates = this%m_coordinates
    constraints = this%m_constraints
    physical = 2 * coordinates
    allocate(position(constraints), rhs(size(y)), gradient(constraints,coordinates), &
        mass(coordinates,coordinates), inverse_gradient(coordinates,constraints), &
        multiplier(constraints,1))
    call this%constraints(x, y(1:coordinates), position, args)
    call this%constraint_jacobian(x, y(1:coordinates), gradient, args)
    multiplier(:,1) = -matmul(gradient, y(coordinates+1:physical)) - this%m_alpha * position
    call ggl_linear_solve(matmul(gradient, transpose(gradient)), multiplier)
    y(physical+constraints+1:) = multiplier(:,1)
    y(physical+1:physical+constraints) = 0.0d0
    call this%compute_rhs(x, y, rhs, args)
    call this%physical_mass(x, y(1:coordinates), mass, args)
    inverse_gradient = transpose(gradient)
    call ggl_linear_solve(mass, inverse_gradient)
    multiplier(:,1) = rhs(physical+constraints+1:)
    call ggl_linear_solve(matmul(gradient, inverse_gradient), multiplier)
    y(physical+1:physical+constraints) = multiplier(:,1)
end subroutine

subroutine ggl_linear_solve(matrix, rhs)
    !! Resolves a square physical or constraint system by pivoted LU, without
    !! forming an inverse. Nonfinite data/invalid shapes cause INVALID_INPUT;
    !! a vanishing or relatively tiny LU pivot causes SINGULAR_MATRIX.
    real(real64), intent(in), dimension(:,:) :: matrix
        !! Nonempty square coefficient matrix, left unchanged.
    real(real64), intent(inout), dimension(:,:) :: rhs
        !! Matching-row matrix of one or more RHS columns, replaced by solutions.
    real(real64), allocatable, dimension(:,:) :: factors
    integer(int32), allocatable, dimension(:) :: pivots
    integer :: equations, equation
    real(real64) :: largest, smallest

    equations = size(matrix,1)
    if (equations < 1 .or. size(matrix,2) /= equations .or. &
        size(rhs,1) /= equations .or. size(rhs,2) < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(matrix)) .or. .not.all(ieee_is_finite(rhs))) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    call lu_factor(matrix, ipvt = pivots, lu = factors)
    if (.not.all(ieee_is_finite(factors))) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    largest = 0.0d0
    smallest = huge(1.0d0)
    do equation = 1, equations
        largest = max(largest, abs(factors(equation,equation)))
        smallest = min(smallest, abs(factors(equation,equation)))
    end do
    if (smallest <= equations * epsilon(1.0d0) * largest) &
        error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    rhs = solve_lu(factors, pivots, rhs)
    if (.not.all(ieee_is_finite(rhs))) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
end subroutine
end module