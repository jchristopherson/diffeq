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

module sparse_chain_model
    !! A planar chain of point masses joined by rigid, massless links.
    !!
    !! Particle 1 is pinned to the origin by the first link, and each later
    !! particle is joined to its predecessor. With P particles there are
    !! N = 2P Cartesian coordinates and M = P holonomic constraints
    !! \(g_k = \|p_k - p_{k-1}\|^2 - L^2 = 0\), where \(p_0 = 0\).
    !!
    !! Each constraint couples only neighbouring particles. As a result, the
    !! mass matrix, the constraint Jacobian, and the augmented GGL Jacobian
    !! are all sparse. [[chain_model]] extends [[ggl_container]] and overrides
    !! the sparse-Jacobian methods. BDF therefore uses its sparse Newton
    !! iteration for the full 2N+2M index-1 system.
    use iso_fortran_env, only : real64
    use diffeq, only : ggl_container
    use linalg, only : csr_matrix, create_csr_matrix
    implicit none
    private
    public :: chain_model
    public :: nlinks, ncoord, nstate
    public :: link_length, particle_mass, gravity
    public :: chain_force, chain_mass, chain_constraints
    public :: chain_gradient, chain_curvature

    integer, parameter :: nlinks = 30
        !! Number of links, particles, and constraints (P = M).
    integer, parameter :: ncoord = 2 * nlinks
        !! Number of Cartesian coordinates N.
    integer, parameter :: nstate = 2 * ncoord + 2 * nlinks
        !! Size of the augmented GGL state [q, v, lambda, mu].
    integer, parameter :: halo = 2
        !! A column owned by particle/constraint j only influences rows owned
        !! by indices j-2 through j+2.
    integer, parameter :: stride = 2 * halo + 1
        !! Minimum owner separation for columns to share a perturbation.
    integer, parameter :: ncolors = 6 * stride
        !! Six column families (q_x, q_y, v_x, v_y, lambda, mu) times stride.
    real(real64), parameter :: link_length = 0.1d0
        !! Length of each link [m].
    real(real64), parameter :: particle_mass = 0.05d0
        !! Mass of each particle [kg].
    real(real64), parameter :: gravity = 9.81d0
        !! Gravitational acceleration [m/s^2].

    type, extends(ggl_container) :: chain_model
        !! A GGL container supplying a sparse augmented Jacobian.
        !!
        !! The Jacobian differentiates the complete GGL right-hand side with
        !! respect to the augmented state. This is required by
        !! [[ggl_container]] when a sparse Jacobian is supplied. Its columns
        !! are computed with central differences of `compute_rhs`.
        !! Structurally orthogonal columns share a perturbation, so each
        !! Jacobian costs 2 * ncolors = 60 RHS evaluations instead of
        !! 2 * nstate = 360.
        integer :: jacobian_calls = 0
            !! The number of sparse-Jacobian evaluations requested by the solver.
    contains
        procedure :: has_sparse_jacobian => chain_has_sparse_jacobian
        procedure :: compute_sparse_jacobian => chain_sparse_jacobian
    end type

contains
! ------------------------------------------------------------------------------
    function chain_has_sparse_jacobian(this) result(rst)
        !! Selects the sparse Newton path in BDF.
        class(chain_model), intent(in) :: this
            !! The chain model.
        logical :: rst
            !! Always true.
        rst = .true.
    end function

! ------------------------------------------------------------------------------
    pure subroutine column_owner(column, owner, family)
        !! Returns the particle or constraint index that owns a state column
        !! and the family (1-6) of that column.
        integer, intent(in) :: column
            !! The augmented-state column index.
        integer, intent(out) :: owner
            !! The owning particle or constraint index.
        integer, intent(out) :: family
            !! 1-2: q_x/q_y, 3-4: v_x/v_y, 5: lambda, 6: mu.

        if (column <= ncoord) then
            owner = (column - 1) / 2 + 1
            family = mod(column - 1, 2) + 1
        else if (column <= 2 * ncoord) then
            owner = (column - ncoord - 1) / 2 + 1
            family = mod(column - ncoord - 1, 2) + 3
        else if (column <= 2 * ncoord + nlinks) then
            owner = column - 2 * ncoord
            family = 5
        else
            owner = column - 2 * ncoord - nlinks
            family = 6
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_sparse_jacobian(this, x, y, jac, args)
        !! Computes the sparse augmented GGL Jacobian with colored central
        !! differences.
        !!
        !! Columns from the same family whose owners differ by at least
        !! `stride` affect disjoint rows. They are perturbed together, and
        !! each nonzero row difference is assigned to the column whose
        !! influence window contains that row.
        class(chain_model), intent(inout) :: this
            !! The chain model.
        real(real64), intent(in) :: x
            !! The independent variable.
        real(real64), intent(in), dimension(:) :: y
            !! The augmented state [q, v, lambda, mu].
        type(csr_matrix), intent(out) :: jac
            !! The sparse Jacobian of the augmented right-hand side.
        class(*), intent(in), optional :: args
            !! Unused user arguments, forwarded to the RHS.

        integer, parameter :: window_rows = (2 * halo + 1) * 6
        real(real64) :: shifted(nstate), plus(nstate), minus(nstate), step(nstate)
        integer :: color_of(nstate), owner_of(nstate)
        integer :: rows(window_rows * nstate), cols(window_rows * nstate)
        real(real64) :: vals(window_rows * nstate)
        integer :: column, color, owner, family, i, k, nnz

        this%jacobian_calls = this%jacobian_calls + 1
        do column = 1, nstate
            call column_owner(column, owner, family)
            owner_of(column) = owner
            color_of(column) = (family - 1) * stride + mod(owner - 1, stride) + 1
            step(column) = 1.0d-6 * max(1.0d0, abs(y(column)))
        end do

        nnz = 0
        do color = 1, ncolors
            shifted = y
            where (color_of == color) shifted = y + step
            call this%compute_rhs(x, shifted, plus, args)
            shifted = y
            where (color_of == color) shifted = y - step
            call this%compute_rhs(x, shifted, minus, args)

            do column = 1, nstate
                if (color_of(column) /= color) cycle
                owner = owner_of(column)
                do i = max(1, owner - halo), min(nlinks, owner + halo)
                    do k = 1, 2
                        call store(2 * (i - 1) + k, column)
                        call store(ncoord + 2 * (i - 1) + k, column)
                    end do
                    call store(2 * ncoord + i, column)
                    call store(2 * ncoord + nlinks + i, column)
                end do
            end do
        end do
        jac = create_csr_matrix(nstate, nstate, rows(1:nnz), cols(1:nnz), vals(1:nnz))

    contains
        subroutine store(row, col)
            integer, intent(in) :: row, col
            real(real64) :: value

            value = (plus(row) - minus(row)) / (2.0d0 * step(col))
            if (value == 0.0d0) return
            nnz = nnz + 1
            rows(nnz) = row
            cols(nnz) = col
            vals(nnz) = value
        end subroutine
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_force(x, y, force, args)
        !! Applied gravitational force; y holds [q, v] of length 2N.
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: force
        class(*), intent(in), optional :: args

        force(1::2) = 0.0d0
        force(2::2) = -particle_mass * gravity
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_mass(x, q, mass, args)
        !! Diagonal physical mass matrix stored in CSR format.
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        type(csr_matrix), intent(out) :: mass
        class(*), intent(in), optional :: args
        integer :: i, index(ncoord)

        index = [(i, i = 1, ncoord)]
        mass = create_csr_matrix(ncoord, ncoord, index, index, &
            spread(particle_mass, 1, ncoord))
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_constraints(x, q, g, args)
        !! Link-length residuals g_k = |p_k - p_{k-1}|^2 - L^2.
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(in), optional :: args
        real(real64) :: d(2)
        integer :: k

        do k = 1, nlinks
            d = link_vector(q, k)
            g(k) = dot_product(d, d) - link_length**2
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_gradient(x, q, gradient, args)
        !! Constraint Jacobian G (M-by-N). Row k contains 2 d_k for particle
        !! k and -2 d_k for particle k-1.
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: gradient
        class(*), intent(in), optional :: args
        real(real64) :: d(2)
        integer :: k

        gradient = 0.0d0
        do k = 1, nlinks
            d = link_vector(q, k)
            gradient(k, 2*k-1:2*k) = 2.0d0 * d
            if (k > 1) gradient(k, 2*k-3:2*k-2) = -2.0d0 * d
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine chain_curvature(x, q, direction, v, curvature, args)
        !! DG(q)[direction] v. The Hessian of g_k is constant, giving
        !! 2 (direction_k - direction_{k-1}) . (v_k - v_{k-1}).
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: q, direction, v
        real(real64), intent(out), dimension(:) :: curvature
        class(*), intent(in), optional :: args
        integer :: k

        do k = 1, nlinks
            curvature(k) = 2.0d0 * dot_product(link_vector(direction, k), &
                link_vector(v, k))
        end do
    end subroutine

! ------------------------------------------------------------------------------
    pure function link_vector(q, k) result(d)
        !! Returns q_k - q_{k-1} for a 2N-vector of particle quantities,
        !! with the anchor (index 0) at zero.
        real(real64), intent(in) :: q(:)
        integer, intent(in) :: k
        real(real64) :: d(2)

        d = q(2*k-1:2*k)
        if (k > 1) d = d - q(2*k-3:2*k-2)
    end function
end module

! ******************************************************************************
program sparse_ggl_example
    !! Simulate a falling 30-link chain with the GGL index reduction and the
    !! sparse BDF solver.
    !!
    !! The chain starts straight, at rest, and at 60 degrees from vertical.
    !! [[ggl_container]] converts the index-3 constrained equations into a
    !! stabilized index-1 DAE. Its sparse mass matrix diag(I, M, 0, 0) is
    !! assembled from the user's CSR physical mass. The overridden
    !! `compute_sparse_jacobian` selects BDF's sparse Newton path.
    !!
    !! The program checks the sparse Jacobian against the dense
    !! finite-difference Jacobian. It also checks that position and velocity
    !! constraint drift stays small.
    use iso_fortran_env, only : real64
    use diffeq, only : bdf
    use linalg, only : csr_matrix, csr_to_dense, nonzero_count
    use fplot_core
    use sparse_chain_model
    implicit none

    ! Parameters
    integer, parameter :: npts = 301
    integer, parameter :: nsnap = 6
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: final_time = 3.0d0
    real(real64), parameter :: release_angle = pi / 3.0d0
    real(real64), parameter :: alpha = 2.0d1, beta = 2.0d1

    ! Model Variables
    type(chain_model) :: model
    type(bdf) :: solver
    type(csr_matrix) :: sparse_jac
    real(real64) :: times(npts), state(nstate), dense_jac(nstate, nstate)
    real(real64) :: g(nlinks), gradient(nlinks, ncoord)
    real(real64) :: position_drift(npts), velocity_drift(npts)
    real(real64) :: jacobian_error, energy0, energy1
    real(real64), allocatable :: sol(:,:)
    integer :: i, j, row

    ! Plot variables
    type(multiplot) :: plt
    type(plot_2d) :: plt_chain, plt_drift
    type(plot_data_2d) :: snapshots(nsnap), pd_position, pd_velocity
    class(plot_axis), pointer :: yaxis
    class(legend), pointer :: lgnd
    real(real64) :: chain_x(nlinks + 1), chain_y(nlinks + 1)
    character(len = 32) :: label

    ! Configure the GGL container. The physical mass is supplied only in sparse
    ! form, and the augmented mass matrix diag(I, M, 0, 0) is built from it.
    call model%initialize(ncoord, nlinks, alpha = alpha, beta = beta)
    model%fcn => chain_force
    model%sparse_mass_matrix => chain_mass
    model%constraints => chain_constraints
    model%constraint_jacobian => chain_gradient
    model%constraint_curvature => chain_curvature
    call model%set_is_mass_matrix_dependent(.false.)

    ! Straight chain at rest. Positions and velocities satisfy the constraints;
    ! initialize_multipliers computes consistent lambda and mu.
    state = 0.0d0
    do j = 1, nlinks
        state(2*j-1) = j * link_length * sin(release_angle)
        state(2*j) = -j * link_length * cos(release_angle)
    end do
    if (.not.model%validate_initial_conditions(0.0d0, state, 1.0d-12)) &
        error stop "Inconsistent initial conditions"
    call model%initialize_multipliers(0.0d0, state)

    ! Check the colored sparse Jacobian against the inherited dense one.
    call model%compute_sparse_jacobian(0.0d0, state, sparse_jac)
    call model%compute_jacobian(0.0d0, state, dense_jac)
    jacobian_error = maxval(abs(csr_to_dense(sparse_jac) - dense_jac)) / &
        max(1.0d0, maxval(abs(dense_jac)))
    print '(A,I0,A,I0,A,I0,A)', "Augmented GGL system: ", nstate, " equations (", &
        ncoord, " coordinates, ", nlinks, " constraints)"
    print '(A,I0,A,I0,A,F5.2,A)', "Jacobian nonzeros: ", nonzero_count(sparse_jac), &
        " of ", nstate**2, " (", 1.0d2 * nonzero_count(sparse_jac) / real(nstate**2, real64), "%)"
    print '(A,ES10.3)', "Relative sparse vs. dense Jacobian difference: ", jacobian_error
    if (jacobian_error > 1.0d-5) error stop "Sparse Jacobian does not match"
    model%jacobian_calls = 0

    ! Integrate
    times = [(final_time * i / (npts - 1.0d0), i = 0, npts - 1)]
    call solver%set_absolute_tolerance(1.0d-7)
    call solver%set_relative_tolerance(1.0d-7)
    call solver%solve(model, times, state)
    sol = solver%get_solution()

    ! Measure constraint drift |g(q)| and |G(q) v| at each output time.
    do row = 1, npts
        call chain_constraints(sol(row,1), sol(row,2:ncoord+1), g)
        call chain_gradient(sol(row,1), sol(row,2:ncoord+1), gradient)
        position_drift(row) = maxval(abs(g))
        velocity_drift(row) = maxval(abs(matmul(gradient, sol(row,ncoord+2:2*ncoord+1))))
    end do
    energy0 = chain_energy(sol(1,2:))
    energy1 = chain_energy(sol(npts,2:))

    print '(A,I0)', "Sparse Jacobian evaluations: ", model%jacobian_calls
    print '(A,ES10.3)', "Maximum position constraint drift: ", maxval(position_drift)
    print '(A,ES10.3)', "Maximum velocity constraint drift: ", maxval(velocity_drift)
    print '(A,ES10.3)', "Relative energy change: ", abs(energy1 - energy0) / abs(energy0)
    if (model%jacobian_calls == 0) error stop "Sparse Jacobian was not used"
    if (maxval(position_drift) > 1.0d-5 .or. maxval(velocity_drift) > 1.0d-4) &
        error stop "Constraint drift exceeded the expected bound"

! ------------------------------------------------------------------------------
    ! Plot chain snapshots and constraint drift.
    call plt%initialize(1, 2, width = 1200)

    call plt_chain%initialize()
    call plt_chain%set_x_axis_title("x [m]")
    call plt_chain%set_y_axis_title("y [m]")
    lgnd => plt_chain%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    do i = 1, nsnap
        row = 1 + (i - 1) * (npts - 1) / (2 * (nsnap - 1))
        chain_x(1) = 0.0d0
        chain_y(1) = 0.0d0
        chain_x(2:) = sol(row, 2:ncoord+1:2)
        chain_y(2:) = sol(row, 3:ncoord+1:2)
        write(label, '("t = ",F4.2)') sol(row,1)
        call snapshots(i)%define_data(chain_x, chain_y)
        call snapshots(i)%set_name(trim(label))
        call snapshots(i)%set_line_width(2.0)
        call snapshots(i)%set_draw_markers(.true.)
        call snapshots(i)%set_marker_style(MARKER_FILLED_CIRCLE)
        call snapshots(i)%set_marker_scaling(0.5)
        call plt_chain%push(snapshots(i))
    end do
    call plt%set(1, 1, plt_chain)

    call plt_drift%initialize()
    call plt_drift%set_x_axis_title("t [s]")
    call plt_drift%set_y_axis_title("Constraint Drift")
    yaxis => plt_drift%get_y_axis()
    call yaxis%set_is_log_scaled(.true.)
    lgnd => plt_drift%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)
    call lgnd%set_horizontal_position(LEGEND_CENTER)
    call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
    ! Exact zeros at t = 0 cannot be displayed on a log axis.
    call pd_position%define_data(sol(:,1), max(position_drift, 1.0d-16))
    call pd_position%set_name("max |g(q)|")
    call pd_position%set_line_width(2.0)
    call pd_velocity%define_data(sol(:,1), max(velocity_drift, 1.0d-16))
    call pd_velocity%set_name("max |G(q)v|")
    call pd_velocity%set_line_width(2.0)
    call plt_drift%push(pd_position)
    call plt_drift%push(pd_velocity)
    call plt%set(1, 2, plt_drift)

    call plt%draw()

contains
    pure function chain_energy(y) result(e)
        !! Kinetic plus gravitational potential energy of [q, v, ...].
        real(real64), intent(in) :: y(:)
        real(real64) :: e

        e = 0.5d0 * particle_mass * sum(y(ncoord+1:2*ncoord)**2) + &
            particle_mass * gravity * sum(y(2:ncoord:2))
    end function
end program
