submodule (diffeq) diffeq_expfixed
    use linalg, only : mtx_mult
contains
! ------------------------------------------------------------------------------
pure module function ef_get_order(this) result(rst)
    class(exponential_fixed_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 3
end function

! ------------------------------------------------------------------------------
module subroutine ef_alloc_workspace(this, neqn, err)
    ! Arguments
    class(exponential_fixed_integrator), intent(inout) :: this
    integer(int32), intent(in) :: neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: flag
    character(len = :), allocatable :: errmsg

    ! Process
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

    if (allocated(this%m_jac2)) then
        if (size(this%m_jac2, 1) /= neqn .or. size(this%m_jac2, 2) /= neqn) then
            deallocate(this%m_jac2)
            allocate(this%m_jac2(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_jac2(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_work)) then
        if (size(this%m_work) /= neqn) then
            deallocate(this%m_work)
            allocate(this%m_work(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_work(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("ef_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine ef_step(this, sys, h, x, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(exponential_fixed_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, j, neqn
    real(real64) :: h2_2, h3_6
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    h2_2 = 0.5d0 * h**2
    h3_6 = (h**2) / 6.0d0
    neqn = size(y)

    ! Input Checking
    if (size(yn) /= neqn) go to 10

    ! Allocate the workspace
    call this%allocate_workspace(neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Evaluate the function
    call sys%fcn(x, y, this%m_work)

    ! Compute the Jacobian (J)
    call sys%compute_jacobian(x, y, this%m_work, this%m_jac, errmgr)
    if (errmgr%has_error_occurred()) return

    ! h**3 / 6 * J**2
    call mtx_mult(.false., .false., h3_6, this%m_jac, this%m_jac, 0.0d0, &
        this%m_jac2)

    ! Compute h * I + (h**2/2) * J + (h**3/6) * J**2 - store in m_jac
    do j = 1, neqn
        do i = 1, neqn
            this%m_jac(i,j) = h2_2 * this%m_jac(i,j) + this%m_jac2(i,j)
            if (i == j) then
                this%m_jac(i,j) = this%m_jac(i,j) + h
            end if
        end do
        yn(j) = y(j)
    end do

    ! Multiply the results by f, and add to y
    call mtx_mult(.false., 1.0d0, this%m_jac, this%m_work, 1.0d0, yn)

    ! End
    return

    ! YN array size error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to have ", neqn, &
        " elements, but was found to have ", size(yn), " elements."
    call errmgr%report_error("ef_step", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule