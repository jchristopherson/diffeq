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
    use linalg, only : identity, csr_matrix, &
        create_empty_csr_matrix, csr_to_dense, dense_to_csr, size, &
        validate_csr_storage
    use linalg_errors, only : LA_SINGULAR_MATRIX_ERROR, LA_CONVERGENCE_ERROR
    use diffeq_sparse_workspace, only : sparse_solver_workspace
    use lapack, only : dgetrf, dgetrs
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
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
            class(*), intent(in), optional :: args
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
            class(*), intent(in), optional :: args
                !! Optional user-owned parameters, forwarded unchanged.
        end subroutine
    end interface

    type :: ggl_linear_workspace
        real(real64), allocatable :: factors(:,:)
        integer, allocatable :: pivots(:)
        type(csr_matrix) :: sparse_mass
        type(sparse_solver_workspace) :: sparse_solver
        real(real64), allocatable :: sparse_solution(:)
    end type

    type :: ggl_rhs_workspace
        real(real64), allocatable :: position(:), curvature(:), qdot(:)
        real(real64), allocatable :: gradient(:,:), mass(:,:), acceleration(:,:)
        type(ggl_linear_workspace) :: linear
    end type

    type :: ggl_initial_workspace
        real(real64), allocatable :: position(:), rhs(:)
        real(real64), allocatable :: gradient(:,:), mass(:,:), inverse_gradient(:,:), multiplier(:,:), constraint_matrix(:,:)
        type(ggl_linear_workspace) :: physical_linear, constraint_linear
    end type

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
        !! The inherited fcn_ad is not supported: differentiating only the
        !! applied force would omit the real-valued constraint assembly.
        !! Assigning fcn_ad causes DIFFEQ_INVALID_OPERATION_ERROR.
        !!
        !! An optional inherited `sparse_mass_matrix` receives q of length N
        !! and writes the physical N-by-N sparse mass. A dense mass callback,
        !! if assigned, takes precedence; otherwise the mass defaults to
        !! identity. Sparse physical mass solves use ILUT-preconditioned GMRES.
        !! The augmented sparse matrix diag(I,M,0,0) is assembled internally.
        !! A sparse `sparse_jacobian` must, like
        !! `jacobian`, differentiate the complete augmented RHS.
        !! RHS scratch storage is retained per OpenMP thread. A container
        !! must not be used concurrently by independent or nested integrations.
        integer, private :: m_coordinates = 0
        integer, private :: m_constraints = 0
        real(real64), private :: m_alpha = 1.0d0
        real(real64), private :: m_beta = 1.0d0
        type(ggl_rhs_workspace), private, allocatable :: m_rhs_work(:)
        type(ggl_initial_workspace), private :: m_initial_work
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
        procedure :: has_sparse_mass_matrix => ggl_has_sparse_mass
        procedure :: get_is_ode_defined => ggl_is_defined
        procedure :: compute_rhs => ggl_rhs
        procedure :: compute_mass_matrix => ggl_mass
        procedure :: compute_sparse_mass_matrix => ggl_sparse_mass
        procedure :: validate_initial_conditions => ggl_validate
        procedure :: initialize_multipliers => ggl_initialize_multipliers
        procedure, private :: check_state => ggl_check_state
        procedure, private :: physical_mass => ggl_physical_mass
    end type

    interface ggl_container
        module procedure :: ggl_constructor_1
    end interface

contains
! ------------------------------------------------------------------------------
function ggl_constructor_1(coordinate_count, constraint_count, alpha, beta) &
    result(rst)
    !! Defines the mechanical dimensions and resets stabilization gains to one
    !! unless supplied. Callback associations and inherited model options are
    !! retained. Invalid counts or gains cause DIFFEQ_INVALID_INPUT_ERROR.
    integer, intent(in) :: coordinate_count
        !! Positive number N of coordinates and physical velocities.
    integer, intent(in) :: constraint_count
        !! Number M of independent constraints, with 1 <= M <= N.
    real(real64), intent(in), optional :: alpha, beta
        !! Nonnegative finite position and velocity decay rates, respectively,
        !! in inverse independent-variable units. Zero disables the respective
        !! damping term without disabling the differentiated constraint.
    type(ggl_container) :: rst
        !! The new [[ggl_container]] object.

    call rst%initialize(coordinate_count, constraint_count, alpha = alpha, &
        beta = beta)
end function

! ------------------------------------------------------------------------------
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
    if (coordinate_count /= this%m_coordinates .or. constraint_count /= this%m_constraints) then
        if (allocated(this%m_rhs_work)) deallocate(this%m_rhs_work)
        this%m_initial_work = ggl_initial_workspace()
    end if
    this%m_coordinates = coordinate_count
    this%m_constraints = constraint_count
    position_rate = 1.0d0
    velocity_rate = 1.0d0
    if (present(alpha)) position_rate = alpha
    if (present(beta)) velocity_rate = beta
    call this%set_stabilization_parameters(position_rate, velocity_rate)
end subroutine

! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
function ggl_state_count(this) result(rst)
    !! Returns 2N+2M, or zero before dimensions have been initialized.
    class(ggl_container), intent(in) :: this
        !! Container whose augmented state size is queried.
    integer :: rst
        !! Required length of the state, initial conditions, and RHS arrays.

    rst = 2 * (this%m_coordinates + this%m_constraints)
end function

! ------------------------------------------------------------------------------
function ggl_has_mass(this) result(rst)
    !! Always returns true: the augmented matrix is singular even when the
    !! physical matrix defaults to identity. Explicit ODE paths are unsuitable.
    class(ggl_container), intent(in) :: this
        !! Container being queried; initialization is not required.
    logical :: rst
        !! True, indicating that a mass-matrix DAE path is required.

    rst = .true.
end function

! ------------------------------------------------------------------------------
function ggl_has_sparse_mass(this) result(rst)
    !! Always returns true: the augmented sparse matrix diag(I,M,0,0) is
    !! assembled internally, even when the physical matrix is identity.
    class(ggl_container), intent(in) :: this
        !! Container being queried; initialization is not required.
    logical :: rst
        !! True, indicating that the sparse mass-matrix path applies.

    rst = .true.
end function

! ------------------------------------------------------------------------------
function ggl_is_defined(this) result(rst)
    !! Checks configuration and all required callback associations. This does
    !! not evaluate callbacks or certify constraint rank or mass regularity.
    class(ggl_container), intent(in) :: this
        !! Container whose model definition is inspected.
    logical :: rst
        !! True only for initialized dimensions and associated force, position,
        !! constraint-Jacobian, and curvature callbacks. Mass is optional.

    if (associated(this%fcn_ad)) error stop DIFFEQ_INVALID_OPERATION_ERROR
    rst = this%m_coordinates > 0 .and. this%m_constraints > 0 .and. &
        associated(this%fcn) .and. associated(this%constraints) .and. &
        associated(this%constraint_jacobian) .and. associated(this%constraint_curvature)
end function

! ------------------------------------------------------------------------------
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

! ------------------------------------------------------------------------------
subroutine ggl_physical_mass(this, x, q, mass, args)
    !! Evaluates the physical mass callback, or supplies identity. A sparse
    !! physical callback is expanded to dense form. This private helper assumes
    !! the caller has checked array dimensions. Singular masses are diagnosed
    !! by the linear solve, not by this evaluation.
    class(ggl_container), intent(inout) :: this
        !! Container holding the optional mechanical mass callback.
    real(real64), intent(in) :: x
        !! Independent-variable value for the mechanical mass callback.
    real(real64), intent(in), dimension(:) :: q
        !! N coordinates, without velocities or multipliers.
    real(real64), intent(out), dimension(:,:) :: mass
        !! Physical N-by-N matrix, not the augmented singular matrix.
    class(*), intent(in), optional :: args
        !! Optional user-owned parameters forwarded to the callback.
    type(csr_matrix) :: sparse

    if (associated(this%mass_matrix)) then
        call this%mass_matrix(x, q, mass, args)
    else if (associated(this%sparse_mass_matrix)) then
        call this%sparse_mass_matrix(x, q, sparse, args)
        if (size(sparse, 1) /= size(q) .or. size(sparse, 2) /= size(q)) &
            error stop DIFFEQ_MATRIX_SIZE_ERROR
        if (validate_csr_storage(sparse) /= 0) error stop DIFFEQ_INVALID_INPUT_ERROR
        mass = csr_to_dense(sparse)
    else
        mass = identity(size(q))
    end if
end subroutine

! ------------------------------------------------------------------------------
subroutine ggl_prepare_rhs_workspace(this, tid)
    class(ggl_container), intent(inout) :: this
    integer, intent(out) :: tid
    integer :: threads, n, m

    tid = omp_get_thread_num() + 1
    threads = omp_get_num_threads()
    n = this%m_coordinates
    m = this%m_constraints
    !$omp critical (ggl_workspace_setup)
    if (allocated(this%m_rhs_work)) then
        if (size(this%m_rhs_work) < threads) deallocate(this%m_rhs_work)
    end if
    if (.not.allocated(this%m_rhs_work)) allocate(this%m_rhs_work(threads))
    if (.not.allocated(this%m_rhs_work(tid)%position)) then
        allocate(this%m_rhs_work(tid)%position(m), this%m_rhs_work(tid)%curvature(m), &
            this%m_rhs_work(tid)%qdot(n), this%m_rhs_work(tid)%gradient(m,n), &
            this%m_rhs_work(tid)%acceleration(n,1))
    end if
    !$omp end critical (ggl_workspace_setup)
end subroutine

! ------------------------------------------------------------------------------
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
    class(*), intent(in), optional :: args
        !! Optional user parameters forwarded unchanged to every callback.
    integer :: coordinates, constraints, physical, tid

    call this%check_state(y)
    if (.not.this%get_is_ode_defined()) error stop DIFFEQ_NULL_POINTER_ERROR
    if (size(dydx) /= size(y)) error stop DIFFEQ_ARRAY_SIZE_ERROR
    coordinates = this%m_coordinates
    constraints = this%m_constraints
    physical = 2 * coordinates
    call ggl_prepare_rhs_workspace(this, tid)
    associate(work => this%m_rhs_work(tid))
    call this%constraints(x, y(1:coordinates), work%position, args)
    call this%constraint_jacobian(x, y(1:coordinates), work%gradient, args)
    work%qdot = y(coordinates+1:physical) + &
        matmul(transpose(work%gradient), y(physical+constraints+1:))
    dydx(1:coordinates) = work%qdot
    call this%fcn(x, y(1:physical), dydx(coordinates+1:physical), args)
    dydx(coordinates+1:physical) = dydx(coordinates+1:physical) - &
        matmul(transpose(work%gradient), y(physical+1:physical+constraints))
    work%acceleration(:,1) = dydx(coordinates+1:physical)
    if (associated(this%sparse_mass_matrix) .and. .not.associated(this%mass_matrix)) then
        call ggl_sparse_physical_solve(this%sparse_mass_matrix, x, y(1:coordinates), &
            work%acceleration, work%linear, args)
    else
        if (.not.allocated(work%mass)) allocate(work%mass(coordinates,coordinates))
        call this%physical_mass(x, y(1:coordinates), work%mass, args)
        call ggl_linear_solve(work%mass, work%acceleration, work%linear)
    end if
    call this%constraint_curvature(x, y(1:coordinates), work%qdot, &
        y(coordinates+1:physical), work%curvature, args)
    dydx(physical+1:physical+constraints) = matmul(work%gradient, work%qdot) + this%m_alpha * work%position
    dydx(physical+constraints+1:) = work%curvature + matmul(work%gradient, work%acceleration(:,1)) + &
        this%m_beta * matmul(work%gradient, y(coordinates+1:physical))
    end associate
end subroutine

! ------------------------------------------------------------------------------
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
    class(*), intent(in), optional :: args
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

! ------------------------------------------------------------------------------
subroutine ggl_sparse_mass(this, x, y, m, args)
    !! Assembles the sparse augmented matrix diag(I,M,0,0) from the physical
    !! mass. Dense and sparse physical callbacks follow the same precedence as
    !! [[ggl_mass]]; without either callback the physical matrix is identity.
    class(ggl_container), intent(inout) :: this
        !! Container defining dimensions and the optional physical sparse mass.
    real(real64), intent(in) :: x
        !! Current independent-variable value.
    real(real64), intent(in), dimension(:) :: y
        !! Augmented state [q,v,lambda,mu] of length 2N+2M.
    type(csr_matrix), intent(out) :: m
        !! Sparse augmented matrix with side length 2N+2M.
    class(*), intent(in), optional :: args
        !! Optional user parameters forwarded to the physical mass callback.
    type(csr_matrix) :: physical
    integer :: coordinates, total, i, nnz
    real(real64), allocatable, dimension(:,:) :: dense

    call this%check_state(y)
    coordinates = this%m_coordinates
    total = size(y)
    if (associated(this%mass_matrix)) then
        allocate(dense(coordinates, coordinates))
        call this%mass_matrix(x, y(1:coordinates), dense, args)
        physical = dense_to_csr(dense)
    else if (associated(this%sparse_mass_matrix)) then
        call this%sparse_mass_matrix(x, y(1:coordinates), physical, args)
    else
        physical = dense_to_csr(identity(coordinates))
    end if
    if (size(physical, 1) /= coordinates .or. size(physical, 2) /= coordinates) &
        error stop DIFFEQ_MATRIX_SIZE_ERROR
    if (validate_csr_storage(physical) /= 0) error stop DIFFEQ_INVALID_INPUT_ERROR
    nnz = size(physical%values)

    m = create_empty_csr_matrix(total, total, coordinates + nnz)
    do i = 1, coordinates
        m%row_indices(i) = i
        m%column_indices(i) = i
        m%values(i) = 1.0d0
        m%row_indices(coordinates+i) = coordinates + physical%row_indices(i)
    end do
    m%row_indices(2*coordinates+1:) = coordinates + nnz + 1
    m%column_indices(coordinates+1:) = coordinates + physical%column_indices
    m%values(coordinates+1:) = physical%values
end subroutine

! ------------------------------------------------------------------------------
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
    class(*), intent(in), optional :: args
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

! ------------------------------------------------------------------------------
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
    class(*), intent(in), optional :: args
        !! Optional user-owned parameters forwarded to model callbacks.
    integer :: coordinates, constraints, physical

    call this%check_state(y)
    if (.not.this%get_is_ode_defined()) error stop DIFFEQ_NULL_POINTER_ERROR
    coordinates = this%m_coordinates
    constraints = this%m_constraints
    physical = 2 * coordinates
    associate(work => this%m_initial_work)
    if (.not.allocated(work%position)) then
        allocate(work%position(constraints), work%rhs(size(y)), &
            work%gradient(constraints,coordinates), work%inverse_gradient(coordinates,constraints), &
            work%multiplier(constraints,1), work%constraint_matrix(constraints,constraints))
    end if
    call this%constraints(x, y(1:coordinates), work%position, args)
    call this%constraint_jacobian(x, y(1:coordinates), work%gradient, args)
    work%multiplier(:,1) = -matmul(work%gradient, y(coordinates+1:physical)) - this%m_alpha * work%position
    work%constraint_matrix = matmul(work%gradient, transpose(work%gradient))
    call ggl_linear_solve(work%constraint_matrix, work%multiplier, work%constraint_linear)
    y(physical+constraints+1:) = work%multiplier(:,1)
    y(physical+1:physical+constraints) = 0.0d0
    call this%compute_rhs(x, y, work%rhs, args)
    work%inverse_gradient = transpose(work%gradient)
    if (associated(this%sparse_mass_matrix) .and. .not.associated(this%mass_matrix)) then
        call ggl_sparse_physical_solve(this%sparse_mass_matrix, x, y(1:coordinates), &
            work%inverse_gradient, work%physical_linear, args)
    else
        if (.not.allocated(work%mass)) allocate(work%mass(coordinates,coordinates))
        call this%physical_mass(x, y(1:coordinates), work%mass, args)
        call ggl_linear_solve(work%mass, work%inverse_gradient, work%physical_linear)
    end if
    work%multiplier(:,1) = work%rhs(physical+constraints+1:)
    work%constraint_matrix = matmul(work%gradient, work%inverse_gradient)
    call ggl_linear_solve(work%constraint_matrix, work%multiplier, work%constraint_linear)
    y(physical+1:physical+constraints) = work%multiplier(:,1)
    end associate
end subroutine

! ------------------------------------------------------------------------------
subroutine ggl_sparse_physical_solve(mass_callback, x, q, rhs, work, args)
    procedure(ode_sparse_mass_matrix) :: mass_callback
    real(real64), intent(in) :: x, q(:)
    real(real64), intent(inout) :: rhs(:,:)
    type(ggl_linear_workspace), intent(inout) :: work
    class(*), intent(in), optional :: args
    integer :: flag, column, n

    n = size(q)
    call mass_callback(x, q, work%sparse_mass, args)
    if (size(work%sparse_mass, 1) /= n .or. size(work%sparse_mass, 2) /= n) &
        error stop DIFFEQ_MATRIX_SIZE_ERROR
    if (.not.allocated(work%sparse_solution)) allocate(work%sparse_solution(n))
    call work%sparse_solver%factor(work%sparse_mass, flag)
    if (flag == LA_SINGULAR_MATRIX_ERROR) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    if (flag /= 0) error stop flag
    do column = 1, size(rhs, 2)
        call work%sparse_solver%solve(work%sparse_mass, rhs(:,column), work%sparse_solution, flag)
        if (flag == LA_SINGULAR_MATRIX_ERROR) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
        if (flag == LA_CONVERGENCE_ERROR) error stop DIFFEQ_CONVERGENCE_ERROR
        if (flag /= 0) error stop flag
        rhs(:,column) = work%sparse_solution
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine ggl_linear_solve(matrix, rhs, work)
    !! Resolves a square physical or constraint system by pivoted LU, without
    !! forming an inverse. Nonfinite data/invalid shapes cause INVALID_INPUT;
    !! a vanishing or relatively tiny LU pivot causes SINGULAR_MATRIX.
    real(real64), intent(in), dimension(:,:) :: matrix
        !! Nonempty square coefficient matrix, left unchanged.
    real(real64), intent(inout), dimension(:,:) :: rhs
        !! Matching-row matrix of one or more RHS columns, replaced by solutions.
    type(ggl_linear_workspace), intent(inout) :: work
    integer :: equations, equation, info
    real(real64) :: largest, smallest

    equations = size(matrix,1)
    if (equations < 1 .or. size(matrix,2) /= equations .or. &
        size(rhs,1) /= equations .or. size(rhs,2) < 1) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(matrix)) .or. .not.all(ieee_is_finite(rhs))) &
        error stop DIFFEQ_INVALID_INPUT_ERROR
    if (allocated(work%pivots)) then
        if (size(work%pivots) /= equations) deallocate(work%pivots, work%factors)
    end if
    if (.not.allocated(work%pivots)) allocate(work%pivots(equations), work%factors(equations,equations))
    work%factors = matrix
    call dgetrf(equations, equations, work%factors, equations, work%pivots, info)
    if (info > 0) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    if (info < 0) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(work%factors))) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    largest = 0.0d0
    smallest = huge(1.0d0)
    do equation = 1, equations
        largest = max(largest, abs(work%factors(equation,equation)))
        smallest = min(smallest, abs(work%factors(equation,equation)))
    end do
    if (smallest <= equations * epsilon(1.0d0) * largest) &
        error stop DIFFEQ_SINGULAR_MATRIX_ERROR
    call dgetrs('N', equations, size(rhs,2), work%factors, equations, work%pivots, rhs, equations, info)
    if (info /= 0) error stop DIFFEQ_INVALID_INPUT_ERROR
    if (.not.all(ieee_is_finite(rhs))) error stop DIFFEQ_SINGULAR_MATRIX_ERROR
end subroutine

! ------------------------------------------------------------------------------
end module