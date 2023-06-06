submodule (diffeq) diffeq_dirk
    ! Model Parameters (ESDIRK4(3)6L[2]SA, Table 16, pg 90)
    ! https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf
    real(real64), parameter :: gamma = 0.25d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = (2.0d0 - sqrt(2.0d0)) / 4.0d0
    real(real64), parameter :: c4 = 5.0d0 / 8.0d0
    real(real64), parameter :: c5 = 26.0d0 / 25.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: a21 = gamma
    real(real64), parameter :: a22 = gamma
    real(real64), parameter :: a32 = (1.0d0 - sqrt(2.0d0)) / 8.0d0
    real(real64), parameter :: a31 = (c3 - a32 - gamma)
    real(real64), parameter :: a33 = gamma
    real(real64), parameter :: a43 = 7.0d0 * (1.0d0 + sqrt(2.0d0)) / 32.0d0
    real(real64), parameter :: a42 = (5.0d0 - 7.0d0 * sqrt(2.0d0)) / 64.0d0
    real(real64), parameter :: a41 = c4 - a42 - a43 - gamma
    real(real64), parameter :: a44 = gamma
    real(real64), parameter :: a54 = &
        166.0d0 * (-97.0d0 + 376.0d0 * sqrt(2.0d0)) / 1.09375d5
    real(real64), parameter :: a53 = &
        (5.06605d5 + 1.32109d5 * sqrt(2.0d0)) / 4.375d5
    real(real64), parameter :: a52 = &
        (-1.3796d4 - 5.4539d4 * sqrt(2.0d0)) / 1.25d5
    real(real64), parameter :: a51 = c5 - a52 - a53 - a54 - gamma
contains
! ------------------------------------------------------------------------------
pure module function dirk_get_order(this) result(rst)
    class(dirk_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

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