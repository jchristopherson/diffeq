submodule (diffeq) diffeq_bsrk32
    real(real64), parameter :: a21 = 0.5d0
    real(real64), parameter :: a31 = 0.0d0
    real(real64), parameter :: a32 = 0.75d0
    real(real64), parameter :: a41 = 2.0d0 / 9.0d0
    real(real64), parameter :: a42 = 1.0d0 / 3.0d0
    real(real64), parameter :: a43 = 4.0d0 / 9.0d0

    real(real64), parameter :: b1 = 2.0d0 / 9.0d0
    real(real64), parameter :: b2 = 1.0d0 / 3.0d0
    real(real64), parameter :: b3 = 4.0d0 / 9.0d0
    real(real64), parameter :: b4 = 0.0d0

    real(real64), parameter :: b1a = 7.0d0 / 24.0d0
    real(real64), parameter :: b2a = 1.0d0 / 4.0d0
    real(real64), parameter :: b3a = 1.0d0 / 3.0d0
    real(real64), parameter :: b4a = 1.0d0 / 8.0d0

    real(real64), parameter :: c1 = 0.0d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = 0.75d0
    real(real64), parameter :: c4 = 1.0d0
contains
! ------------------------------------------------------------------------------
module subroutine bsrk32_define_model(this)
    ! Arguments
    class(bsrk32_integrator), intent(inout) :: this

    ! Process
    if (this%m_modelDefined) return

    ! A
    this%a = 0.0d0

    this%a(2,1) = a21

    this%a(3,1) = a31
    this%a(3,2) = a32

    this%a(4,1) = a41
    this%a(4,2) = a42
    this%a(4,3) = a43

    ! B
    this%b(1) = b1
    this%b(2) = b2
    this%b(3) = b3
    this%b(4) = b4

    ! C
    this%c(1) = c1
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4

    ! E
    this%e(1) = b1a - b1
    this%e(2) = b2a - b2
    this%e(3) = b3a - b3
    this%e(4) = b4a - b4
end subroutine

! ------------------------------------------------------------------------------
pure module function bsrk32_is_fsal(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_stage_count(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 4
end function

! ------------------------------------------------------------------------------
pure module function bsrk32_get_order(this) result(rst)
    class(bsrk32_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 3
end function

! ------------------------------------------------------------------------------
module subroutine bsrk32_set_up_interp(this, x, xn, y, yn, k)
    ! Arguments
    class(bsrk32_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn
    real(real64), intent(in), dimension(:,:) :: k

    ! Parameters
    integer(int32), parameter :: n = 3

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h

    ! Initialization
    neqn = size(y)
    h = this%get_step_size()

    ! Memory ALlocation
    if (allocated(this%m_bsrk23work)) then
        if (size(this%m_bsrk23work, 1) /= neqn .or. &
            size(this%m_bsrk23work, 2) /= n) &
        then
            deallocate(this%m_bsrk23work)
            allocate(this%m_bsrk23work(neqn, n))
        end if
    else
        allocate(this%m_bsrk23work(neqn, n))
    end if

    ! Construct the coefficient arrays
    this%m_bsrk23work(:,1) = -(y - yn + k(:,1) * h) / h**2
    this%m_bsrk23work(:,2) = k(:,1) - 2.0d0 * x * this%m_bsrk23work(:,1)
    this%m_bsrk23work(:,3) = y - this%m_bsrk23work(:,1) * x**2 - &
        this%m_bsrk23work(:,2) * x
end subroutine

! ------------------------------------------------------------------------------
module subroutine bsrk32_interp(this, xprev, xnew, x, y, err)
    ! Arguments
    class(bsrk32_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: neqn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(this%m_bsrk23work, 1)

    ! Input Check
    if (size(y) /= neqn) go to 10

    ! Process
    y = this%m_bsrk23work(:,1) * x**2 + this%m_bsrk23work(:,2) * x + &
        this%m_bsrk23work(:,3)

    ! End
    return

    ! Y is not sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to be of length ", &
        neqn, ", but was found to be of length ", size(y), "."
    call errmgr%report_error("bsrk32_interp", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule