submodule (diffeq) diffeq_rkfixed
contains
! ------------------------------------------------------------------------------
! n = method order
! neqn = # of ODE's to solve
module subroutine rkf_alloc_workspace(this, n, neqn, err)
    ! Arguments
    class(rk_fixed_integrator), intent(inout) :: this
    integer(int32), intent(in) :: n, neqn
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: flag, mw, nw
    character(len = :), allocatable :: errmsg

    ! Process
    mw = neqn
    nw = n + 1  ! +1 for 1 additional NEQN-element workspace array
    if (allocated(this%m_work)) then
        if (size(this%m_work, 1) /= mw .or. size(this%m_work, 2) /= nw) then
            deallocate(this%m_work)
            allocate(this%m_work(mw, nw), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_work(mw, nw), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("rkf_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rkf_step(this, sys, h, x, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(rk_fixed_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, j, n, neqn, m, n1
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(y)
    n = this%get_order()
    n1 = n + 1  ! index of additional workspace array #1

    ! Input Checking
    if (size(yn) /= neqn) go to 10

    ! Allocate the workspace
    call this%allocate_workspace(n, neqn, errmgr)
    if (errmgr%has_error_occurred()) return

    ! As this is an explicit routine, the Butcher tableau is lower triangular.
    call sys%fcn(x, y, this%m_work(:,1))
    do i = 2, n
        this%m_work(:,n1) = 0.0d0
        do j = 1, i - 1 ! only reference the sub-diagonal components
            this%m_work(:,n1) = this%m_work(:,n1) + &
                this%get_method_factor(i,j) * this%m_work(:,j)
        end do

        call sys%fcn( &
            x + h * this%get_position_factor(i), &
            y + h * this%m_work(:,n1), &
            this%m_work(:,i) &   ! output
        )
    end do

    ! Compute the next solution estimate
    this%m_work(:,n1) = 0.0d0
    do i = 1, n
        this%m_work(:,n1) = this%m_work(:,n1) + &
            this%get_quadrature_weight(i) * this%m_work(:,i)
    end do
    yn = y + h * this%m_work(:,n1)

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