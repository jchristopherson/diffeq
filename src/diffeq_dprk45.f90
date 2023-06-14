submodule (diffeq) diffeq_dprk45
    ! Dormand-Prince 4th/5th Order Model Coefficients
    real(real64), parameter :: a21 = 1.0d0 / 5.0d0
    real(real64), parameter :: a31 = 3.0d0 / 40.0d0
    real(real64), parameter :: a32 = 9.0d0 / 40.0d0
    real(real64), parameter :: a41 = 44.0d0 / 45.0d0
    real(real64), parameter :: a42 = -56.0d0 / 15.0d0
    real(real64), parameter :: a43 = 32.0d0 / 9.0d0
    real(real64), parameter :: a51 = 1.9372d4 / 6.561d3
    real(real64), parameter :: a52 = -2.536d4 / 2.187d3
    real(real64), parameter :: a53 = 6.4448d4 / 6.561d3
    real(real64), parameter :: a54 = -2.12d2 / 7.29d2
    real(real64), parameter :: a61 = 9.017d3 / 3.168d3
    real(real64), parameter :: a62 = -3.55d2 / 33.0d0
    real(real64), parameter :: a63 = 4.6732d4 / 5.247d3
    real(real64), parameter :: a64 = 49.0d0 / 1.76d2
    real(real64), parameter :: a65 = -5.103d3 / 1.8656d4
    real(real64), parameter :: a71 = 35.0d0 / 3.84d2
    real(real64), parameter :: a72 = 0.0d0
    real(real64), parameter :: a73 = 5.0d2 / 1.113d3
    real(real64), parameter :: a74 = 1.25d2 / 1.92d2
    real(real64), parameter :: a75 = -2.187d3 / 6.784d3
    real(real64), parameter :: a76 = 11.0d0 / 84.0d0
    
    real(real64), parameter :: e1 = -71.0d0 / 5.76d4
    real(real64), parameter :: e2 = 0.0d0
    real(real64), parameter :: e3 = 71.0d0 / 1.6695d4
    real(real64), parameter :: e4 = -71.0d0 / 1.92d3
    real(real64), parameter :: e5 = 1.7253d4 / 3.392d5
    real(real64), parameter :: e6 = -22.0d0 / 5.25d2
    real(real64), parameter :: e7 = 1.0d0 / 4.0d1

    real(real64), parameter :: c2 = 1.0d0 / 5.0d0
    real(real64), parameter :: c3 = 3.0d0 / 1.0d1
    real(real64), parameter :: c4 = 4.0d0 / 5.0d0
    real(real64), parameter :: c5 = 8.0d0 / 9.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: c7 = 1.0d0

    ! Interpolation Parameters
    real(real64), parameter :: d1 = -1.2715105075d10 / 1.1282082432d10
    real(real64), parameter :: d3 = 8.74874797d10 / 3.2700410799d10
    real(real64), parameter :: d4 = -1.0690763975d10 / 1.880347072d9
    real(real64), parameter :: d5 = 7.01980252875d11 / 1.99316789632d11
    real(real64), parameter :: d6 = -1.453857185d9 / 8.22651844d8
    real(real64), parameter :: d7 = 6.9997945d7 / 2.9380423d7
contains
! ------------------------------------------------------------------------------
module subroutine dprk45_define_model(this)
    ! Arguments
    class(dprk45_integrator), intent(inout) :: this

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

    this%a(5,1) = a51
    this%a(5,2) = a52
    this%a(5,3) = a53
    this%a(5,4) = a54

    this%a(6,1) = a61
    this%a(6,2) = a62
    this%a(6,3) = a63
    this%a(6,4) = a64
    this%a(6,5) = a65
    
    this%a(7,1) = a71
    this%a(7,2) = a72
    this%a(7,3) = a73
    this%a(7,4) = a74
    this%a(7,5) = a75
    this%a(7,6) = a76

    ! B
    this%b(1) = a71
    this%b(2) = a72
    this%b(3) = a73
    this%b(4) = a74
    this%b(5) = a75
    this%b(6) = a76
    this%b(7) = 0.0d0

    ! C
    this%c(1) = 0.0d0
    this%c(2) = c2
    this%c(3) = c3
    this%c(4) = c4
    this%c(5) = c5
    this%c(6) = c6
    this%c(7) = c7

    ! E
    this%e(1) = e1
    this%e(2) = e2
    this%e(3) = e3
    this%e(4) = e4
    this%e(5) = e5
    this%e(6) = e6
    this%e(7) = e7

    ! Update definition status
    this%m_modelDefined = .true.
end subroutine

! ------------------------------------------------------------------------------
pure module function dprk45_is_fsal(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    logical :: rst
    rst = .true.
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_stage_count(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 7
end function

! ------------------------------------------------------------------------------
pure module function dprk45_get_order(this) result(rst)
    class(dprk45_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 5
end function

! ------------------------------------------------------------------------------
module subroutine dprk45_set_up_interp(this, x, xn, y, yn, k)
    ! Arguments
    class(dprk45_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn
    real(real64), intent(in), dimension(:,:) :: k

    ! Parameters
    integer(int32), parameter :: n = 5

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h, ydiff, bspl

    ! Intialization
    neqn = size(y)
    h = this%get_step_size()

    ! Memory Allocation
    if (allocated(this%m_dprk45work)) then
        if (size(this%m_dprk45work, 1) /= neqn .or. &
            size(this%m_dprk45work, 2) /= n) &
        then
            deallocate(this%m_dprk45work)
            allocate(this%m_dprk45work(neqn, n))
        end if
    else
        allocate(this%m_dprk45work(neqn, n))
    end if

    ! Construct the coefficient arrays
    do i = 1, neqn
        ydiff = yn(i) - y(i)
        bspl = h * k(i,1) - ydiff

        this%m_dprk45work(i,1) = y(i)
        this%m_dprk45work(i,2) = ydiff
        this%m_dprk45work(i,3) = bspl
        this%m_dprk45work(i,4) = ydiff - h * k(i,7) - bspl
        this%m_dprk45work(i,5) = h * (d1 * k(i,1) + d3 * k(i,3) + &
            d4 * k(i,4) + d5 * k(i,5) + d6 * k(i,6) + d7 * k(i,7))
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine dprk45_interp(this, xprev, yprev, xnew, x, y, err)
    ! Arguments
    class(dprk45_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(in), dimension(:) :: yprev
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: i, neqn
    real(real64) :: h, s, s1
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(this%m_dprk45work, 1)
    h = xnew - xprev
    s = (x - xprev) / h
    s1 = 1.0d0 - s

    ! Input Check
    if (size(y) /= neqn) go to 10

    ! Process
    y = this%m_dprk45work(:,1) + s * ( &
        this%m_dprk45work(:,2) + s1 * ( &
        this%m_dprk45work(:,3) + s * ( &
        this%m_dprk45work(:,4) + s1 * this%m_dprk45work(:,5) &
    )))

    ! End
    return

    ! Y is not sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to be of length ", &
        neqn, ", but was found to be of length ", size(y), "."
    call errmgr%report_error("dprk45_interp", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule