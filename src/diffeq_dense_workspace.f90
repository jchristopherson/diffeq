! This file is part of diffeq.
!
! diffeq is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! diffeq is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with diffeq. If not, see <https://www.gnu.org/licenses/>.

module diffeq_dense_workspace
    !! Reusable dense linear-system storage for the implicit integrators.
    !!
    !! A system is factored with pivoted LU unless QR is explicitly requested.
    !! LU is rejected if its reciprocal one-norm condition estimate is no
    !! larger than machine epsilon. Each LU solution is also checked against
    !! a true residual scaled by \(100\epsilon(\|A\|_\infty\|x\|_\infty+
    !! \|b\|_\infty)\). A rejected LU uses LAPACK DGELSY's
    !! rank-revealing complete orthogonal factorization and minimum-norm
    !! least-squares solution, with epsilon as the relative rank threshold.
    !! No factors are reused after a new matrix is set. A workspace must not
    !! be used concurrently by multiple solves.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use lapack, only : dgetrf, dgetrs, dgelsy
    implicit none
    private
    public :: dense_solver_workspace

    interface
        subroutine dgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
            !! LAPACK condition estimator; the bundled lapack module does not expose it.
            import :: real64
            character, intent(in) :: norm
            integer, intent(in) :: n, lda
            real(real64), intent(in) :: a(lda,*), anorm
            real(real64), intent(out) :: rcond, work(*)
            integer, intent(out) :: iwork(*), info
        end subroutine
    end interface

    type :: dense_solver_workspace
        !! Scratch and factorizations for one square system at a time.
        !! Buffers are retained until the matrix dimension changes.
        private
        real(real64), allocatable :: matrix(:,:), factors(:,:), qr(:,:)
            !! Original matrix, LU factors, and disposable QR input.
        real(real64), allocatable :: rhs(:), work(:), condition_work(:)
            !! QR right-hand side, LAPACK work, and condition-estimation work.
        integer, allocatable :: pivots(:), qr_pivots(:), condition_iwork(:)
            !! LU pivots, QR column pivots, and condition-estimation work.
        logical :: use_qr = .true.
            !! True if the current matrix requires the rank-revealing path.
        logical :: ready = .false.
            !! True only after a successful factorization request.
    contains
        procedure :: factor => dense_factor
            !! Stores a new system and prepares LU where appropriate.
        procedure :: solve => dense_solve
            !! Solves a right-hand side using the currently stored system.
    end type

contains
subroutine dense_factor(this, matrix, info, qr_only)
    !! Initializes or resizes scratch storage and prepares the supplied matrix.
    !! `info` is zero on success, including when QR fallback is selected.
    class(dense_solver_workspace), intent(inout) :: this
        !! Workspace whose previous matrix and factors are invalidated.
    real(real64), intent(in) :: matrix(:,:)
        !! Nonempty, finite square coefficient matrix.
    integer, intent(out) :: info
        !! Zero on success, otherwise a LAPACK or invalid-input status.
    logical, intent(in), optional :: qr_only
        !! If true, retain QR least-squares behavior even for nonsingular matrices.
    integer :: n, j, rank, lwork
    real(real64) :: anorm, rcond, query(1)

    this%ready = .false.
    n = size(matrix,1)
    info = -1
    if (n < 1 .or. size(matrix,2) /= n) return
    if (.not.all(ieee_is_finite(matrix))) return
    if (allocated(this%pivots)) then
        if (size(this%pivots) /= n) deallocate(this%matrix, this%factors, this%qr, &
            this%rhs, this%work, this%condition_work, this%pivots, &
            this%qr_pivots, this%condition_iwork)
    end if
    if (.not.allocated(this%pivots)) then
        allocate(this%matrix(n,n), this%factors(n,n), this%qr(n,n), this%rhs(n), &
            this%pivots(n), this%qr_pivots(n), this%condition_work(4*n), &
            this%condition_iwork(n))
        this%qr = 0.0d0
        this%rhs = 0.0d0
        this%qr_pivots = 0
        call dgelsy(n, n, 1, this%qr, n, this%rhs, n, this%qr_pivots, &
            epsilon(1.0d0), rank, query, -1, info)
        if (info /= 0 .or. .not.ieee_is_finite(query(1))) error stop "DGELSY workspace query failed"
        lwork = max(1, ceiling(query(1)))
        allocate(this%work(lwork))
    end if
    this%matrix = matrix
    this%use_qr = .false.
    if (present(qr_only)) this%use_qr = qr_only
    info = 0
    if (this%use_qr) then
        this%ready = .true.
        return
    end if
    this%factors = matrix
    call dgetrf(n, n, this%factors, n, this%pivots, info)
    if (info < 0) return
    if (info > 0) then
        this%use_qr = .true.
        info = 0
        this%ready = .true.
        return
    end if
    if (.not.all(ieee_is_finite(this%factors))) then
        this%use_qr = .true.
        this%ready = .true.
        return
    end if
    anorm = 0.0d0
    do j = 1, n
        anorm = max(anorm, sum(abs(matrix(:,j))))
    end do
    if (.not.ieee_is_finite(anorm)) then
        this%use_qr = .true.
        this%ready = .true.
        return
    end if
    call dgecon('1', n, this%factors, n, anorm, rcond, &
        this%condition_work, this%condition_iwork, info)
    if (info /= 0) return
    this%use_qr = .not.ieee_is_finite(rcond)
    if (ieee_is_finite(rcond)) this%use_qr = rcond <= epsilon(1.0d0)
    this%ready = .true.
end subroutine

subroutine dense_solve(this, b, x, info)
    !! Solves the current system without allocating per-call scratch.
    !! LU solutions must satisfy a scaled backward-error check; otherwise
    !! DGELSY computes a rank-revealing minimum-norm least-squares solution.
    class(dense_solver_workspace), intent(inout) :: this
        !! Workspace previously prepared by `factor`.
    real(real64), intent(in) :: b(:)
        !! Finite right-hand side with one entry per matrix row.
    real(real64), intent(out) :: x(:)
        !! Solution, with one entry per matrix column.
    integer, intent(out) :: info
        !! Zero on success, otherwise a LAPACK or invalid-input status.
    integer :: n, i, rank
    real(real64) :: error_norm, scale
    logical :: fallback

    n = size(b)
    info = -1
    if (.not.this%ready) return
    if (.not.allocated(this%pivots)) return
    if (size(this%pivots) /= n .or. size(x) /= n .or. .not.all(ieee_is_finite(b))) return
    fallback = this%use_qr
    if (.not.fallback) then
        x = b
        call dgetrs('N', n, 1, this%factors, n, this%pivots, x, n, info)
        if (info /= 0) return
        fallback = .not.all(ieee_is_finite(x))
        if (.not.fallback) then
            error_norm = 0.0d0
            scale = maxval(abs(b))
            do i = 1, n
                error_norm = max(error_norm, abs(dot_product(this%matrix(i,:), x) - b(i)))
                scale = max(scale, sum(abs(this%matrix(i,:))) * maxval(abs(x)) + maxval(abs(b)))
            end do
            fallback = .not.ieee_is_finite(error_norm) .or. .not.ieee_is_finite(scale)
            if (.not.fallback) fallback = error_norm > 100.0d0 * epsilon(1.0d0) * scale
        end if
    end if
    if (fallback) then
        this%use_qr = .true.
        this%qr = this%matrix
        this%rhs = b
        this%qr_pivots = 0
        call dgelsy(n, n, 1, this%qr, n, this%rhs, n, this%qr_pivots, &
            epsilon(1.0d0), rank, this%work, size(this%work), info)
        if (info /= 0) return
        x = this%rhs
        if (.not.all(ieee_is_finite(x))) info = -1
    end if
end subroutine
end module
