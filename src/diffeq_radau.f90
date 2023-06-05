submodule (diffeq) diffeq_radau
    use linalg

    ! Model Constants
    real(real64), parameter :: sq6 = sqrt(6.d0)
    real(real64), parameter :: c1 = (4.d0 - sq6) / 10.d0
    real(real64), parameter :: c2 = (4.d0 + sq6) / 10.d0
    real(real64), parameter :: c1m1 = c1 - 1.d0
    real(real64), parameter :: c2m1 = c2 - 1.d0
    real(real64), parameter :: c1mc2 = c1 - c2
    real(real64), parameter :: dd1 =  - (13.d0 + 7.d0*sq6) / 3.d0
    real(real64), parameter :: dd2 = ( - 13.d0 + 7.d0*sq6) / 3.d0
    real(real64), parameter :: dd3 =  - 1.d0 / 3.d0
    real(real64), parameter :: u1 = 30.0d0 / &
        (6.0d0 + 81.0d0**(1.0d0 / 3.0d0) - 9.0d0**(1.0d0 / 3.0d0))
    real(real64), parameter :: alph = (12.0d0 - 81.0d0**(1.0d0 / 3.0d0) + &
        9.0d0**(1.0d0 / 3.0d0)) / 60.0d0
    real(real64), parameter :: beta = (81.0d0**(1.0d0 / 3.0d0) + &
        9.0d0**(1.0d0 / 3.0d0)) * sqrt(3.0d0) / 60.0d0
    real(real64), parameter :: cno = alph**2 + beta**2
    real(real64), parameter :: t11 = 9.1232394870892942792d-02
    real(real64), parameter :: t12 =  - 0.14125529502095420843d0
    real(real64), parameter :: t13 =  - 3.0029194105147424492d-02
    real(real64), parameter :: t21 = 0.24171793270710701896d0
    real(real64), parameter :: t22 = 0.20412935229379993199d0
    real(real64), parameter :: t23 = 0.38294211275726193779d0
    real(real64), parameter :: t31 = 0.96604818261509293619d0
    real(real64), parameter :: ti11 = 4.3255798900631553510d0
    real(real64), parameter :: ti12 = 0.33919925181580986954d0
    real(real64), parameter :: ti13 = 0.54177053993587487119d0
    real(real64), parameter :: ti21 =  - 4.1787185915519047273d0
    real(real64), parameter :: ti22 =  - 0.32768282076106238708d0
    real(real64), parameter :: ti23 = 0.47662355450055045196d0
    real(real64), parameter :: ti31 =  - 0.50287263494578687595d0
    real(real64), parameter :: ti32 = 2.5719269498556054292d0
    real(real64), parameter :: ti33 =  - 0.59603920482822492497d0
contains
! ------------------------------------------------------------------------------
module subroutine rad_alloc_workspace(this, neqn, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
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

    ! Call the base process
    call vsi_alloc_workspace(this, neqn, errmgr)
    if (errmgr%has_error_occurred()) return

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

    if (allocated(this%m_mass)) then
        if (size(this%m_mass, 1) /= neqn .or. &
            size(this%m_mass, 2) /= neqn) &
        then
            deallocate(this%m_mass)
            allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_mass(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_e1)) then
        if (size(this%m_e1, 1) /= neqn .or. size(this%m_e1, 2) /= neqn) then
            deallocate(this%m_e1)
            allocate(this%m_e1(neqn, neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_e1(neqn, neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_e2)) then
        if (size(this%m_e2, 1) /= neqn .or. size(this%m_e2, 2) /= neqn) then
            deallocate(this%m_e2)
            allocate(this%m_e2(neqn, neqn), stat = flag, &
                source = (0.0d0, 0.0d0))
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_e2(neqn, neqn), stat = flag, source = (0.0d0, 0.0d0))
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_ie1)) then
        if (size(this%m_ie1) /= neqn) then
            deallocate(this%m_ie1)
            allocate(this%m_ie1(neqn), stat = flag, source = 0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_ie1(neqn), stat = flag, source = 0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_ie2)) then
        if (size(this%m_ie2) /= neqn) then
            deallocate(this%m_ie2)
            allocate(this%m_ie2(neqn), stat = flag, source = 0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_ie2(neqn), stat = flag, source = 0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_cont)) then
        if (size(this%m_cont) /= 4 * neqn) then
            deallocate(this%m_cont)
            allocate(this%m_cont(4 * neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_cont(4 * neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_z1)) then
        if (size(this%m_z1) /= neqn) then
            deallocate(this%m_z1)
            allocate(this%m_z1(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_z1(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_z2)) then
        if (size(this%m_z2) /= neqn) then
            deallocate(this%m_z2)
            allocate(this%m_z2(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_z2(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_z3)) then
        if (size(this%m_z3) /= neqn) then
            deallocate(this%m_z3)
            allocate(this%m_z3(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_z3(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_zc)) then
        if (size(this%m_zc) /= neqn) then
            deallocate(this%m_zc)
            allocate(this%m_zc(neqn), stat = flag, source = (0.0d0, 0.0d0))
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_zc(neqn), stat = flag, source = (0.0d0, 0.0d0))
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_f1)) then
        if (size(this%m_f1) /= neqn) then
            deallocate(this%m_f1)
            allocate(this%m_f1(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_f1(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_f2)) then
        if (size(this%m_f2) /= neqn) then
            deallocate(this%m_f2)
            allocate(this%m_f2(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_f2(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_f3)) then
        if (size(this%m_f3) /= neqn) then
            deallocate(this%m_f3)
            allocate(this%m_f3(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_f3(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_scale)) then
        if (size(this%m_scale) /= neqn) then
            deallocate(this%m_scale)
            allocate(this%m_scale(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_scale(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if

    if (allocated(this%m_dydx)) then
        if (size(this%m_dydx) /= neqn) then
            deallocate(this%m_dydx)
            allocate(this%m_dydx(neqn), stat = flag, source = 0.0d0)
            if (flag /= 0) go to 10
        end if
    else
        allocate(this%m_dydx(neqn), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
    end if


    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call err%report_error("rad_alloc_workspace", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_build_e1(this, fac, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in) :: fac
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, j, neqn
    
    ! Initialization
    neqn = size(this%m_jac, 1)

    ! Process
    if (this%m_usemass) then
        this%m_e1 = this%m_mass * fac - this%m_jac
    else
        do j = 1, neqn
            do i = 1, neqn
                this%m_e1(i,j) = -this%m_jac(i,j)
            end do
            this%m_e1(j,j) = this%m_e1(j,j) + fac
        end do
    end if

    ! Compute the LU factorization of E1
    call lu_factor(this%m_e1, this%m_ie1, err)
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_build_e2(this, alphan, betan, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in) :: alphan, betan
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, j, neqn
    
    ! Initialization
    neqn = size(this%m_jac, 1)

    ! Process
    if (this%m_usemass) then
        this%m_e2 = cmplx( &
            this%m_mass * alphan - this%m_jac, &
            this%m_mass * betan, &
            real64 &
        )
    else
        do j = 1, neqn
            do i = 1, neqn
                if (i == j) then
                    this%m_e2(i,j) = cmplx( &
                        alphan - this%m_jac(i,j), &
                        betan, &
                        real64 &
                    )
                else
                    this%m_e2(i,j) = cmplx(-this%m_jac(i,j), 0.0d0, real64)
                end if
            end do
        end do
    end if

    ! Compute the LU factorization of E2
    call lu_factor(this%m_e2, this%m_ie2, err)
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_use_default_newton(this) result(rst)
    class(radau_integrator), intent(in) :: this
    logical :: rst
    rst = this%m_usedefaultstart
end function

! --------------------
module subroutine rad_set_use_default_newton(this, x)
    class(radau_integrator), intent(inout) :: this
    logical, intent(in) :: x
    this%m_usedefaultstart = x
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_index_2(this) result(rst)
    class(radau_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_index2
end function

! --------------------
module subroutine rad_set_index_2(this, x)
    class(radau_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_index2 = x
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_index_3(this) result(rst)
    class(radau_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_index3
end function

! --------------------
module subroutine rad_set_index_3(this, x)
    class(radau_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_index3 = x
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_max_newton_iter(this) result(rst)
    class(radau_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = this%m_maxnewton
end function

! --------------------
module subroutine rad_set_max_newton_iter(this, x)
    class(radau_integrator), intent(inout) :: this
    integer(int32), intent(in) :: x
    this%m_maxnewton = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_set_up_newton(this, first, neqn, h)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    logical, intent(in) :: first
    integer(int32), intent(in) :: neqn
    real(real64), intent(in) :: h

    ! Local Variables
    integer(int32) :: i, n2, n3
    real(real64) :: c3q, c1q, c2q, hold, ak1, ak2, ak3

    ! Initialization
    hold = this%get_previous_step_size()
    n2 = 2 * neqn
    n3 = 3 * neqn

    ! Process
    if (first .or. this%get_use_default_newton_start()) then
        do i = 1, neqn
            this%m_z1(i) = 0.0d0
            this%m_z2(i) = 0.0d0
            this%m_z3(i) = 0.0d0
            this%m_f1(i) = 0.0d0
            this%m_f2(i) = 0.0d0
            this%m_f3(i) = 0.0d0
        end do
    else
        c3q = h / hold
        c1q = c1 * c3q
        c2q = c2 * c3q
        do i = 1, neqn
            ak1 = this%m_cont(i+neqn)
            ak2 = this%m_cont(i+n2)
            ak3 = this%m_cont(i+n3)
            z1i = c1q * (ak1 + (c1q - c2m1) * (ak2 + (c1q - c1m1) * ak3))
            z2i = c2q * (ak1 + (c2q - c2m1) * (ak2 + (c2q - c1m1) * ak3))
            z3i = c3q * (ak1 + (c3q - c2m1) * (ak2 + (c3q - c1m1) * ak3))
            this%m_z1(i) = z1i
            this%m_z2(i) = z2i
            this%m_z3(i) = z3i
            this%m_f1(i) = ti11 * z1i + ti12 * z2i + ti13 * z3i
            this%m_f2(i) = ti21 * z1i + ti22 * z2i + ti23 * z3i
            this%m_f3(i) = ti31 * z1i + ti32 * z2i + ti33 * z3i
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_solve_linear_system(this, fac, alphan, betan)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in) :: fac, alphan, betan

    ! Local Variables
    integer(int32) :: i, j, neqn
    real(real64) :: s1, s2, s3, bb

    ! Initialization
    neqn = size(this%m_f1)

    ! Process
    if (this%m_usemass) then
        do i = 1, neqn
            s1 = 0.0d0
            s2 = 0.0d0
            s3 = 0.0d0
            do j = 1, neqn
                bb = this%m_mass(i,j)
                s1 = s1 - bb * this%m_f1(j)
                s2 = s2 - bb * this%m_f2(j)
                s3 = s3 - bb * this%m_f3(j)
            end do
            this%m_z1(i) = this%m_z1(i) + s1 * fac
            this%m_z2(i) = this%m_z2(i) + s2 * alphan - s3 * betan
            this%m_z3(i) = this%m_z3(i) + s3 * alpha + s2 * betan
            this%m_zc(i) = cmplx(this%m_z2(i), this%m_z3(i), real64)
        end do
    else
        do i = 1, neqn
            s2 = -this%m_f2(i)
            s3 = -this%m_f3(i)
            this%m_z1(i) = this%m_z1(i) - this%m_f1(i) * fac
            this%m_z2(i) = this%m_z2(i) + s2 * alphan - s3 * betan
            this%m_z3(i) = this%m_z3(i) + s3 * alphan + s2 * betan
            this%m_zc(i) = cmplx(this%m_z2(i), this%m_z3(i), real64)
        end do
    end if

    ! Solve the two linear systems
    call solve_lu(this%m_e1, this%m_ie1, this%m_z1)
    call solve_lu(this%m_e2, this%m_ie2, this%m_zc)
    do i = 1, neqn
        this%m_z2(i) = real(this%m_zc(i))
        this%m_z3(i) = aimag(this%m_zc(i))
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_estrad(this, sys, h, x, y, first, reject, en)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    logical, intent(in) :: first, reject
    real(real64), intent(out), dimension(:) :: en

    ! Local Variables
    integer(int32) :: i, j, neqn
    real(real64) :: hee1, hee2, hee3, sum

    ! Initialization
    neqn = size(y)
    hee1 = dd1 / h
    hee2 = dd2 / h
    hee3 = dd3 / h

    ! Process
    if (this%m_usemass) then
        this%m_f1 = hee1 * this%m_z1 + hee2 + this%m_z2 + hee3 * this%m_z3
        do i = 1, neqn
            sum = 0.0d0
            do j = 1, neqn
                sum = sum + this%m_mass(i,j) * this%m_f1(j)
            end do
            this%m_f2(i) = sum
            this%m_cont(i) = sum + y(i)
        end do
    else
        do i = 1, neqn
            this%m_f2(i) = hee1 * this%m_z1(i) + hee2 * this%m_z2(i) + &
                hee3 * this%m_z3(i)
            this%m_cont(i) = this%m_f2(i) + y(i)
        end do
    end if
    call solve_lu(this%m_e1, this%m_ie1, this%m_cont(1:neqn))

    err = 0.0d0
    do i = 1, neqn
        en(i) = this%m_cont(i) / this%m_scale(i)
        err = err + en(i)**2
    end do
    err = max(sqrt(err / neqn), 1.0d-10)
    if (err < 1.0d0) return
    if (first .or. reject) then
        this%m_cont(1:neqn) = y + this%m_cont(1:neqn)
        call sys%ode(x, this%m_cont(1:neqn), this%m_f1)
        this%m_cont(1:neqn) = this%m_f1 + this%m_f2
        call solve_lu(this%m_e1, this%m_ie1, this%m_cont(1:neqn))
        err = 0.0d0
        do i = 1, neqn
            err = err + (this%m_cont(i) / this%m_scale(i))**2
        end do
        err = max(sqrt(err / neqn), 1.0d-10)
    end if
end subroutine

! ------------------------------------------------------------------------------
module recursive subroutine rad_newton_iteration(this, sys, h, x, y, dydx, &
    niter, reject, first, last, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(inout) :: h
    real(real64), intent(in) :: x
    real(real64), intent(in), dimension(:) :: y, dydx
    integer(int32), intent(inout) :: niter
    logical, intent(inout) :: reject, first, last
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: neqn, maxiter
    real(real64) :: fac1, alphn, betan
    logical :: converged
    character(len = :), allocatable :: errmsg

    ! Initialization
    neqn = size(y)
    maxiter = this%get_max_newton_iteration()

    ! Update the Jacobian?
    if (this%m_updatejac) then
        call sys%compute_jacobian(x, y, dydx, this%m_jac, err)
        if (err%has_error_occurred()) return
        this%m_updatejac = .false.
    end if

    ! Compute and factor matrices E1 & E2
    fac1 = u1 / h
    alphn = alph / h
    betan = beta / h
    call this%build_e1(fac1, err)
    if (err%has_error_occurred()) return

    call this%build_e2(alphn, betan, err)
    if (err%has_error_occurred()) return

    !     if (0.1d0 * abs(h) <= abs(x) * uround) then
    !         ! TO DO: GO TO 177
    !     end if
    !     if (index2) then
    !         ! TO DO
    !     end if
    !     if (index3) then
    !         ! TO DO
    !     end if

    ! Set up the starting values for the Newton iteration
    call this%set_up_newton_iteration(first, neqn, h)

    ! Newton Iteration
    do
        ! Check against iteration limits.  This operation is necessary here
        ! to avoid any potential issues with a recursive call
        if (niter > maxiter) go to 10

        ! Compute the right-hand-side
        this%m_cont(1:neqn) = y + this%m_z1
        call sys%ode(x + c1 * h, this%m_cont(1:neqn), this%m_z1)

        this%m_cont(1:neqn) = y + this%m_z2
        call sys%ode(x + c2 * h, this%m_cont(1:neqn), this%m_z2)

        this%m_cont(1:neqn) = y + this%m_z3
        call sys%ode(x + h, this%m_cont(1:neqn), this%m_z3)

        ! Solve the linear system
        call this%solve_linear_systems(fac1, alphn, betan)

        ! Check for convergence
        converged = this%newton_converged(niter, h, reject, last)
        if (reject) call this%newton_iteration(sys, h, x, y, dydx, niter, &
            reject, first, last, err)
        if (err%has_error_occurred()) return
        if (converged) exit
    end do

    ! End
    return

    ! Newton Iteration Count Exceeded
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Newton iteration count exceeded at x = ", x, "."
    call err%report_error("rad_newton_iteration", trim(errmsg), &
        DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)
    return

    ! Formatting
100 format(A, E10.3, A)
end subroutine

! ------------------------------------------------------------------------------
module function rad_is_newton_converged(this, niter, h, reject, last) &
    result(rst)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    integer(int32), intent(inout) :: niter
    real(real64), intent(inout) :: h
    logical, intent(inout) :: reject, last
    logical :: rst

    ! Local Variables
    integer(int32) :: i, neqn, n3, maxiter
    real(real64) :: dyno, dynold, denom, thq, thqold, theta, faccon, fnewt, &
        dyth, qnewt, hhfac, uround, f1i, f2i, f3i

    ! Initialization
    neqn = size(this%m_scale)
    n3 = 3 * neqn
    maxiter = this%get_max_newton_iteration()
    faccon = 1.0d0
    fnewt = this%get_newton_tolerance()
    uround = epsilon(uround)
    
    ! Process
    rst = .false.
    niter = niter + 1
    dyno = 0.0d0
    do i = 1, neqn
        denom = this%m_scale(i)
        dyno = dyno + (this%m_z1(i) / denom)**2 + (this%m_z2(i) / denom)**2 + &
            (this%m_z3(i) / denom)**2
    end do
    dyno = sqrt(dyno / n3)

    if (niter > 1 .and. niter < maxiter) then
        thq = dyno / dynold
        if (niter == 2) then
            theta = thq
        else
            theta = sqrt(thq * thqold)
        end if
        thqold = thq
        if (theta < 0.99d0) then
            faccon = theta / (1.0d0 - theta)
            dyth = faccon * dyno * theta**(maxiter - 1 - niter) / fnewt
            if (dyth >= 1.0d0) then
                qnewt = max(1.0d-4, min(2.0d1, dyth))
                hhfac = 0.8d0 * qnewt**(-1.0d0 / (4.0d0 + maxiter - 1 - niter))
                h = hhfac * h
                reject = .true.
                last = .false.
            end if
        else
            ! TO DO: GO TO 78
        end if
    end if
    dynold = max(dyno, uround)
    do i = 1, neqn
        f1i = this%m_f1(i) + this%m_z1(i)
        f2i = this%m_f2(i) + this%m_z2(i)
        f3i = this%m_f3(i) + this%m_z3(i)
        this%m_f1(i) = f1i
        this%m_f2(i) = f2i
        this%m_f3(i) = f3i
        this%m_z1(i) = t11 * f1i + t12 * f2i + t13 * f3i
        this%m_z2(i) = t21 * f1i + t22 * t2i + t23 * f3i
        this%m_z3(i) = t31 * f1i + f2i
    end do
    if (faccon * dyno > fnewt) then
        ! No convergence, 
        rst = .false.
    else
        ! Newton iterations have converged
        rst = .true.
    end if
end function

! ------------------------------------------------------------------------------
module subroutine rad_attempt_step(this, sys, h, x, y, yn, en, xprev, yprev, &
    fprev, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: h, x
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn, en
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: neqn, niter
    logical :: updatejac, reject, first, last
    real(real64) :: hn
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    neqn = size(y)
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    niter = 0
    reject = .false.
    first = .false.
    last = .false.
    hn = h

    ! Evaluate and solve the nonlinear system using a Newton iteration scheme
    call this%newton_iteration(sys, hn, x, y, this%m_dydx, niter, &
        reject, first, last, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute an estimate to the error
    call this%estimate_error(sys, hn, x, y, first, reject, en)

    ! Update y
    yn = y + this%m_z3

    ! End
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_step(this, sys, x, xmax, y, yn, xprev, yprev, fprev, err)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in) :: x, xmax
    real(real64), intent(in), dimension(:) :: y
    real(real64), intent(out), dimension(:) :: yn
    real(real64), intent(in), optional, dimension(:) :: xprev
    real(real64), intent(in), optional, dimension(:,:) :: yprev
    real(real64), intent(inout), optional, dimension(:,:) :: fprev
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: neqn, nstep, maxstep
    real(real64) :: h, enorm
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
    maxstep = this%get_max_per_step_iteration_count()
    this%m_scale = this%m_atol + abs(y) * this%m_rtol
    h = this%get_next_step_size()

    ! Input Checking
    if (size(yn) /= neqn) go to 10

    ! Compute the mass matrix, if needed
    ! TO DO: Update
    this%m_usemass = .false.

    ! Update the derivative estimate
    call sys%ode(x, y, this%m_dydx)

    ! Process
    nstep = 0
    do
        ! Attempt a step
        call this%attempt_step(sys, h, x, y, yn, this%m_ework, &
            xprev, yprev, fprev, errmgr)

        ! Compute the error
        enorm = max(sqrt(sum(this%m_ework**2) / neqn), 1.0d-10)

        ! Compute a new step size
        h = this%compute_next_step_size(h, enorm, this%m_enormPrev)

        ! Check to see if the step size is too small
        if (abs(h) < abs(this%get_min_step_size())) go to 20

        ! Do we need to limit the step size to not overstep a limiting x value
        if (this%get_respect_x_max() .and. &
            abs(x + h) > abs(xmax)) &
        then
            h = xmax - x
        end if

        ! Is the step successful (enorm is normalized to the error tolerances
        ! such that a value less or equal to 1 is successful; else, keep going)
        if (enorm <= 1.0d0) exit

        ! Update the iteration counter
        nstep = nstep + 1
        if (nstep > this%get_max_per_step_iteration_count()) go to 30
    end do

    ! Store error values from this step
    this%m_enormPrev = enorm

    ! Store the updated step size
    call this%set_next_step_size(h)
    
    ! Perform any actions needed on a successful step
    call this%on_successful_step(x, x + this%get_step_size(), y, yn)

    ! End
    return

    ! YN array size error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to have ", neqn, &
        " elements, but was found to have ", size(yn), " elements."
    call errmgr%report_error("rad_step", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Step size too small error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "A step size of ", h, " was encountered at x = ", x, &
        ", and is too small to continue integration."
    call errmgr%report_error("rad_step", trim(errmsg), &
        DIFFEQ_STEP_SIZE_TOO_SMALL_ERROR)
    return

    ! Iteration count exceeded error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 102) "The allowable iteration count was exceeded at x = ", &
        x, "."
    call errmgr%report_error("rad_step", trim(errmsg), &
        DIFFEQ_ITERATION_COUNT_EXCEEDED_ERROR)
    return

! Formatting
100 format(A, I0, A, I0, A)
101 format(A, E10.3, A, E10.3, A)
102 format(A, E10.3, A)
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_newton_tol(this) result(rst)
    class(radau_integrator), intent(in) :: this
    real(real64) :: rst
    rst = this%m_newtontol
end function

! --------------------
module subroutine rad_set_newton_tol(this, x)
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in) :: x
    this%m_newtontol = x
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_on_successful_step(this, x, xn, y, yn)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in) :: x, xn
    real(real64), intent(in), dimension(:) :: y, yn

    ! Local Variables
    integer(int32) :: neqn

    ! Initialization
    neqn = size(y)

    ! Process
    this%m_cont(1:neqn) = yn
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_set_up_interp(this, y)
    ! Arguments
    class(radau_integrator), intent(inout) :: this
    real(real64), intent(in), dimension(:) :: y

    ! Local Variables
    integer(int32) :: i, neqn, n2, n3
    real(real64) :: z1i, z2i, ak, acont3

    ! Initialization
    neqn = size(y)
    n2 = 2 * neqn
    n3 = 3 * neqn

    ! Process
    do i = 1, neqn
        z2i = this%m_z2(i)
        z1i = this%m_z1(i)
        this%m_cont(i + neqn) = (z2i - this%m_z3(i)) / c2m1
        ak = (z1i - z2i) / c1mc2
        acont3 = z1i / c1
        acont3 = (ak - acont3) / c2
        this%m_cont(i + n2) = (ak - this%m_cont(i + neqn)) / c1m1
        this%m_cont(i + n3) = this%m_cont(i + n2) - acont3
    end do
end subroutine

! ------------------------------------------------------------------------------
module subroutine rad_interp(this, xprev, xnew, x, y, err)
    ! Arguments
    class(radau_integrator), intent(in) :: this
    real(real64), intent(in) :: xprev, xnew, x
    real(real64), intent(out), dimension(:) :: y
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: n, np1, n2, n2p1, n3, n3p1
    real(real64) :: s, h
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(this%m_f1)
    np1 = n + 1
    n2 = 2 * n
    n2p1 = n2 + 1
    n3 = 3 * n
    n3p1 = n3 + 1
    h = xnew - xprev
    s = (x - xprev) / h

    ! Input Check
    if (size(y) /= neqn) go to 10

    ! Process
    y = this%m_cont(1:n) + s * (this%m_cont(np1:n2) + &
        (s - c2m1) * (this%m_cont(n2p1:n3) + (s - c1m1) * this%m_cont(n3p1:)))

    ! End
    return

    ! Y is not sized correctly
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The output array was expected to be of length ", &
        neqn, ", but was found to be of length ", size(y), "."
    call errmgr%report_error("rad_interp", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
pure module function rad_get_order(this) result(rst)
    class(radau_integrator), intent(in) :: this
    integer(int32) :: rst
    rst = 5
end function

! ------------------------------------------------------------------------------
end submodule