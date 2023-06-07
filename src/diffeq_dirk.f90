submodule (diffeq) diffeq_dirk
    
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