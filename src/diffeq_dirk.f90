submodule (diffeq) diffeq_dirk
    use linalg
contains
! ------------------------------------------------------------------------------
module subroutine dirk_alloc_workspace(this, neqn, err)
    ! Arguments
    class(dirk_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: flag
    character(len = :), allocatable :: errmsg
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Allocations
    if (allocated(this%m_jac)) then
        if (size(this%m_jac, 1) /= neqn .or. size(this%m_jac, 2) /= neqn) then
            deallocate(this%m_jac)
            allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_jac(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_mass)) then
        if (size(this%m_mass, 1) /= neqn .or. size(this%m_mass, 2) /= neqn) then
            deallocate(this%m_mass)
            allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_mtx)) then
        if (size(this%m_mtx, 1) /= neqn .or. size(this%m_mtx, 2) /= neqn) then
            deallocate(this%m_mtx)
            allocate(this%m_mtx(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_mtx(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_pvt)) then
        if (size(this%m_pvt) /= neqn) then
            deallocate(this%m_pvt)
            allocate(this%m_pvt(neqn), stat = flag, source = 0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_pvt(neqn), stat = flag, source = 0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("dirk_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine dirk_build_factored_matrix(this, sys, h, x, y, err)
    ! Arguments
    class(dirk_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    logical :: change
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    change = .false.

    ! Compute the Jacobian, if an update is needed.
    if (.not.this%get_is_jacobian_current() .or. this%m_firstStep) then
        ! Recompute the Jacobian
        call sys%compute_jacobian(x, y, this%m_jac, errmgr)
        if (errmgr%has_error_occurred()) return
        call this%set_is_jacobian_current(.true.)
        change = .true.
    end if

    ! Compute the mass matrix
    if (associated(sys%mass_matrix) .and. &
        (sys%get_is_mass_matrix_dependent() .or. this%m_firstStep)) &
    then
        call sys%mass_matrix(x, y, this%m_mass)
        change = .true.
    end if

    ! Compute the system matrix
    if (associated(sys%mass_matrix)) then
        call this%build_newton_matrix(h, this%m_jac, this%m_mtx, this%m_mass)
    else
        call this%build_newton_matrix(h, this%m_jac, this%m_mtx)
    end if

    ! Compute the LU factorization of the system matrix if anything has changed
    if (change) then
        call lu_factor(this%m_mtx, this%m_pvt, errmgr)
        if (errmgr%has_error_occurred()) return
    end if
end subroutine

! ------------------------------------------------------------------------------
module subroutine dirk_build_matrix(this, h, jac, x, m, err)
    ! Arguments
    class(dirk_integrator), intent(in) :: this
    real(real64), intent(in) :: h
    real(real64), intent(in), dimension(:,:) :: jac
    real(real64), intent(out), dimension(:,:) :: x
    real(real64), intent(in), dimension(:,:), optional :: m
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, j, neqn
    real(real64) :: fac
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg

    ! Initialization
    neqn = size(this%m_jac, 1)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Checking
    if (size(jac, 1) /= neqn .or. size(jac, 2) /= neqn) go to 10
    if (size(x, 1) /= neqn .or. size(x, 2) /= neqn) go to 20
    if (present(m)) then
        if (size(m, 1) /= neqn .or. size(m, 2) /= neqn) go to 30
    end if

    ! Process
    fac = 1.0d0 / (this%get_method_factor(1,1) * h)
    if (present(m)) then
        x = m * fac - jac
    else
        do j = 1, neqn
            do i = 1, neqn
                if (i == j) then
                    x(i,j) = fac - jac(i,j)
                else
                    x(i,j) = -jac(i,j)
                end if
            end do
        end do
    end if

    ! End
    return

    ! Jacobian matrix is not sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The Jacobian matrix was expected to be ", neqn, &
        "-by-", neqn, ", but was found to be ", size(jac, 1), "-by-", &
        size(jac, 2), "."
    call errmgr%report_error("dirk_build_matrix", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! Output matrix is not sized correctly
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output matrix was expected to be ", neqn, &
        "-by-", neqn, ", but was found to be ", size(x, 1), "-by-", &
        size(x, 2), "."
    call errmgr%report_error("dirk_build_matrix", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! Mass matrix is not sized correctly
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The mass matrix was expected to be ", neqn, &
        "-by-", neqn, ", but was found to be ", size(m, 1), "-by-", &
        size(m, 2), "."
    call errmgr%report_error("dirk_build_matrix", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
! module subroutine dirk_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
!     fprev, err)
!     ! Arguments
!     class(dirk_integrator), intent(inout) :: this
!     class(ode_container), intent(inout) :: sys
!     real(real64), intent(in) :: h, x
!     real(real64), intent(in), dimension(:) :: y
!     real(real64), intent(out), dimension(:) :: yn, en
!     real(real64), intent(in), optional, dimension(:) :: xprev
!     real(real64), intent(in), optional, dimension(:,:) :: yprev
!     real(real64), intent(inout), optional, dimension(:,:) :: fprev
!     class(errors), intent(inout), optional, target :: err

!     ! Local Variables
!     class(errors), pointer :: errmgr
!     type(errors), target :: deferr
    
!     ! Initialization
!     if (present(err)) then
!         errmgr => err
!     else
!         errmgr => deferr
!     end if

!     !
! end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule