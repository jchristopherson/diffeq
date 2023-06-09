submodule (diffeq) diffeq_ode_container
contains
! ------------------------------------------------------------------------------
pure module function oc_get_fd_step(this) result(rst)
    class(ode_container), intent(in) :: this
    real(real64) :: rst
    rst = this%m_fdStep
end function

! --------------------
module subroutine oc_set_fd_step(this, x)
    class(ode_container), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_fdStep = x
end subroutine

! ------------------------------------------------------------------------------
pure module function oc_get_is_mass_dependent(this) result(rst)
    class(ode_container), intent(in) :: this
    logical :: rst
    rst = this%m_massDependent
end function

! --------------------
module subroutine oc_set_is_mass_dependent(this, x)
    class(ode_container), intent(inout) :: this
    logical :: x
    this%m_massDependent = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine oc_jacobian(this, x, y, jac, err)
    ! Arguments
    class(ode_container), intent(inout) :: this
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:,:) :: jac
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, ndof
    real(real64) :: h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    ndof = size(y)
    h = this%get_finite_difference_step()

    ! Input Checking
    if (size(jac, 1) /= ndof .or. size(jac, 2) /= ndof) go to 10

    ! Use a user-defined routine, and then be done
    if (associated(this%jacobian)) then
        call this%jacobian(x, y, jac)
        return
    end if

    ! Allocate workspace.  No action is taken if the proper workspace is
    ! already allocated.
    call this%allocate_workspace(ndof, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Finite Difference Approximation
    ! J(i,j) = df(i) / dy(j)
    this%m_jwork(1:ndof) = y
    call this%ode(x, y, this%m_jwork(ndof+1:))
    do i = 1, ndof
        this%m_jwork(i) = this%m_jwork(i) + h
        call this%ode(x, this%m_jwork(1:ndof), jac(:,i))
        jac(:,i) = (jac(:,i) - this%m_jwork(ndof+1:)) / h
        this%m_jwork(i) = y(i)
    end do

    ! End
    return

    ! Jacobian Size Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output matrix must be (", n, "-by-", n, &
        "), but was found to be (", size(jac, 1), "-by-", size(jac, 2), ")."
    call errmgr%report_error("oc_jacobian", trim(errmsg), &
        DIFFEQ_MATRIX_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine oc_alloc_workspace(this, ndof, err)
    ! Arguments
    class(ode_container), intent(inout) :: this
    integer(int32), intent(in) :: ndof
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: flag
    character(len = :), allocatable :: errmsg

    ! Jacobian Workspace Allocation
    if (allocated(this%m_jwork)) then
        if (size(this%m_jwork) /= 2 * ndof) then
            deallocate(this%m_jwork)
            allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_jwork(2 * ndof), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("oc_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine oc_ode_fcn(this, x, y, dydx)
    ! Arguments
    class(ode_container), intent(in) :: this
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: dydx

    ! Process
    if (associated(this%fcn)) call this%fcn(x, y, dydx)
end subroutine

! ------------------------------------------------------------------------------
end submodule