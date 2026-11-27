module diffeq_sparse_workspace
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use linalg, only : csr_matrix, size
    use linalg_errors, only : LA_NO_ERROR, LA_INVALID_INPUT_ERROR, &
        LA_MATRIX_FORMAT_ERROR, LA_MATRIX_NOT_SQUARE_ERROR, &
        LA_NOT_INITIALIZED_ERROR, LA_SINGULAR_MATRIX_ERROR, &
        LA_CONVERGENCE_ERROR, LA_INVALID_OPERATION_ERROR
    use sparskit, only : ilut, pgmres
    implicit none
    private
    public :: sparse_solver_workspace

    type :: sparse_solver_workspace
        private
        integer, allocatable :: ju(:), jlu(:), jw(:), seen(:)
        real(real64), allocatable :: alu(:), w(:), basis(:,:), residual(:)
    contains
        procedure :: factor => sparse_factor
        procedure :: solve => sparse_solve
    end type

contains
subroutine sparse_factor(this, a, flag)
    class(sparse_solver_workspace), intent(inout) :: this
    type(csr_matrix), intent(in) :: a
    integer, intent(out) :: flag
    integer :: n, row, k, lfil, capacity, ierr, stored

    flag = LA_NOT_INITIALIZED_ERROR
    if (.not.allocated(a%row_indices) .or. .not.allocated(a%column_indices) .or. &
        .not.allocated(a%values)) return
    n = size(a, 1)
    flag = LA_INVALID_INPUT_ERROR
    if (n < 1 .or. .not.all(ieee_is_finite(a%values))) return
    flag = LA_MATRIX_NOT_SQUARE_ERROR
    if (size(a, 2) /= n) return
    stored = size(a%values)
    flag = LA_MATRIX_FORMAT_ERROR
    if (size(a%column_indices) /= stored) return
    if (a%row_indices(1) /= 1 .or. a%row_indices(n+1) /= stored+1) return
    if (any(a%row_indices < 1) .or. any(a%row_indices > stored+1)) return
    if (any(a%row_indices(2:) < a%row_indices(:n))) return
    if (any(a%column_indices < 1) .or. any(a%column_indices > n)) return

    if (allocated(this%ju)) then
        if (size(this%ju) /= n) deallocate(this%ju, this%jw, this%seen, &
            this%w, this%basis, this%residual)
    end if
    if (.not.allocated(this%ju)) then
        allocate(this%ju(n), this%jw(2*n), this%seen(n), this%w(n+1), &
            this%basis(n,min(n,50)+1), this%residual(n))
    end if
    this%seen = 0
    lfil = 1
    do row = 1, n
        do k = a%row_indices(row), a%row_indices(row+1)-1
            if (this%seen(a%column_indices(k)) == row) return
            this%seen(a%column_indices(k)) = row
        end do
        if (all(a%values(a%row_indices(row):a%row_indices(row+1)-1) == 0.0d0)) then
            flag = LA_SINGULAR_MATRIX_ERROR
            return
        end if
        lfil = max(lfil, a%row_indices(row+1)-a%row_indices(row))
    end do
    capacity = max(n+1, lfil*n+1, stored)
    if (allocated(this%alu)) then
        if (size(this%alu) < capacity) deallocate(this%alu, this%jlu)
    end if
    if (.not.allocated(this%alu)) allocate(this%alu(capacity), this%jlu(capacity))
    do
        capacity = size(this%alu)
        ! The linked SPARSKIT build does not promise reentrant local storage.
        !$omp critical (diffeq_sparskit)
        call ilut(n, a%values, a%column_indices, a%row_indices, lfil, &
            epsilon(1.0d0), this%alu, this%jlu, this%ju, capacity, &
            this%w, this%jw, ierr)
        !$omp end critical (diffeq_sparskit)
        if (ierr == 0) exit
        if (ierr == -2 .or. ierr == -3) then
            capacity = capacity + max(2*n, capacity/2)
            deallocate(this%alu, this%jlu)
            allocate(this%alu(capacity), this%jlu(capacity))
        else if (ierr == -4) then
            lfil = n
        else if (ierr > 0 .or. ierr == -5) then
            flag = LA_SINGULAR_MATRIX_ERROR
            return
        else if (ierr == -1) then
            flag = LA_MATRIX_FORMAT_ERROR
            return
        else
            flag = LA_INVALID_OPERATION_ERROR
            return
        end if
    end do
    flag = LA_SINGULAR_MATRIX_ERROR
    if (.not.all(ieee_is_finite(this%alu(1:n)))) return
    if (any(this%alu(1:n) == 0.0d0)) return
    if (.not.all(ieee_is_finite(this%alu(n+2:this%jlu(n+1)-1)))) return
    flag = LA_NO_ERROR
end subroutine

subroutine sparse_solve(this, a, b, x, flag)
    class(sparse_solver_workspace), intent(inout) :: this
    type(csr_matrix), intent(in) :: a
    real(real64), intent(in) :: b(:)
    real(real64), intent(out) :: x(:)
    integer, intent(out) :: flag
    integer :: n, row, k, ierr
    real(real64) :: tol

    n = size(b)
    flag = LA_INVALID_INPUT_ERROR
    if (size(x) /= n .or. .not.all(ieee_is_finite(b))) return
    x = 0.0d0
    flag = LA_NO_ERROR
    if (all(b == 0.0d0)) return
    tol = sqrt(epsilon(1.0d0))
    this%residual = b
    !$omp critical (diffeq_sparskit)
    call pgmres(n, min(n,50), this%residual, x, this%basis, tol, 200, 0, &
        a%values, a%column_indices, a%row_indices, this%alu, this%jlu, this%ju, ierr)
    !$omp end critical (diffeq_sparskit)
    this%residual = b
    do row = 1, n
        do k = a%row_indices(row), a%row_indices(row+1)-1
            this%residual(row) = this%residual(row) - a%values(k) * x(a%column_indices(k))
        end do
    end do
    if ((ierr /= 0 .and. ierr /= -1) .or. .not.all(ieee_is_finite(x)) .or. &
        .not.all(ieee_is_finite(this%residual))) then
        flag = LA_CONVERGENCE_ERROR
    else if (norm2(this%residual) > tol * norm2(b)) then
        flag = LA_CONVERGENCE_ERROR
    end if
end subroutine
end module
