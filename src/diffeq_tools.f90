!
module diffeq_tools
    use iso_fortran_env
    use diffeq
    use ferror
    use nonlin_core
    use nonlin_solve
    use linalg
    implicit none
    private
    public :: find_equilibrium_points
    public :: is_stable_equilibrium

contains
! ------------------------------------------------------------------------------
!> @brief Finds equilibrium points for each differential equation.
!!
!! @par Remarks
!! This routine utilizes an iterative solver to attempt to find the roots of 
!! the system of autonomous differential equations y'(x) = f(y) = 0.  The solver
!! utilized by default is a Quasi-Newton type that is capable of returning only
!! a single root per equation; therefore, if multiple roots exist for an 
!! equation, the starting guess must be altered accordingly in an effort to 
!! search for the root.  With that said, the user can choose to use a different
!! solver; however, this routine is still only set-up to accomodate a single
!! root per equation.
!!
!! @param[in,out] sys The system of differential equations to analyze.
!! @param[in] xi An NEQN-element array containing an initial guess at the root
!!  locations.
!! @param[in,out] solver An optional argument, that if supplied, lets the user
!!  specify the solver being used.  The default solver is a Quasi-Newton type
!!  that will use an update strategy regarding the Jacobian to minimize the
!!  number of Jacobian evaluations necessary.
!! @param[in,out] err An optional errors-based object that if provided 
!!  can be used to retrieve information relating to any errors 
!!  encountered during execution. If not provided, a default 
!!  implementation of the errors class is used internally to provide 
!!  error handling.
!!
!! @return The routine returns an NEQN-by-2 matrix with the first column 
!!  containing the root values, and the second column containing the values
!!  of the differential equations at the root values.
function find_equilibrium_points(sys, xi, solver, err) result(rst)
    ! Arguments
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: xi
    class(equation_solver), intent(inout), optional, target :: solver
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable, dimension(:,:) :: rst

    ! Local Variables
    integer(int32) :: neqn, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    class(equation_solver), pointer :: solverptr
    type(quasi_newton_solver), target :: defsolver
    type(vecfcn_helper) :: container
    procedure(vecfcn), pointer :: eqns
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(solver)) then
        solverptr => solver
    else
        solverptr => defsolver
    end if
    neqn = size(xi)
    eqns => eqn2solve
    call container%set_fcn(eqns, neqn, neqn)

    ! Memory Allocation
    allocate(rst(neqn, 2), stat = flag, source = 0.0d0)
    if (flag /= 0) go to 10

    ! Solve the system of equations
    rst(:,1) = xi   ! <- initial guess
    call solverptr%solve(container, rst(:,1), rst(:,2), err = errmgr)

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("find_equilibrium_points", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Formatting
100 format(A, I0, A)

contains
    ! The vector-valued function to solve
    subroutine eqn2solve(x_, f_)
        ! Arguments
        real(real64), intent(in), dimension(:) :: x_
        real(real64), intent(out), dimension(:) :: f_

        ! As this should be an autonomous system (i.e. the only appearance of
        ! the independent variable is in the derivative), the independent 
        ! variable can be set to any number.  A value of 0 is chosen here for
        ! convenience.
        call sys%ode(0.0d0, x_, f_)
    end subroutine
end function

! ------------------------------------------------------------------------------
!> @brief Tests to see if the supplied equilibrium point is stable. 
!!
!! @par Remarks
!! An equilibrium point is considered stable if all of it's eigenvalues have
!! a negative-valued real component.  If any of the eigenvalues have a 
!! positive-valued real component the equilibrium point is unstable.  If all
!! eigenvalues have a positive-valued real component the equilibrium point is
!! referred to as a source, or an unstable focus in the event of complex-valued
!! eigenvalues.  If only some of the eigenvalues have positive-valued real
!! components, the equilibrium point is a saddle point.
!!
!! @param[in,out] sys The system of differential equations to analyze.
!! @param[in] xeq An NEQN-element array containing the equilibrium point to 
!!  analyze.
!! @param[out] lambda An optional NEQN-element array that, if supplied, can be
!!  used to retrieve the eigenvalues of the system Jacobian about the 
!!  equilibrium point @p xeq.
!! @param[in,out] err An optional errors-based object that if provided 
!!  can be used to retrieve information relating to any errors 
!!  encountered during execution. If not provided, a default 
!!  implementation of the errors class is used internally to provide 
!!  error handling.
!!
!! @return Returns true if the equilibrium point is stable; else, returns false.
function is_stable_equilibrium(sys, xeq, lambda, err) result(rst)
    ! Arguments
    class(ode_container), intent(inout) :: sys
    real(real64), intent(in), dimension(:) :: xeq
    complex(real64), intent(out), dimension(:), target, optional :: lambda
    class(errors), intent(inout), optional, target :: err
    logical :: rst

    ! Local Variables
    integer(int32) :: i, neqn, flag
    real(real64), allocatable, dimension(:,:) :: jac
    complex(real64), allocatable, dimension(:), target :: defvals
    complex(real64), pointer, dimension(:) :: vals
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    rst = .false.
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    neqn = size(xeq)
    if (present(lambda)) then
        if (size(lambda) /= neqn) go to 20
        vals(1:neqn) => lambda
    else
        allocate(defvals(neqn), stat = flag)
        if (flag /= 0) go to 10
        vals(1:neqn) => defvals
    end if

    ! Memory Allocation
    allocate(jac(neqn, neqn), stat = flag, source = 0.0d0)
    if (flag /= 0) go to 10

    ! Determine the Jacobian, and it's eigenvalues
    call sys%compute_jacobian(0.0d0, xeq, jac, errmgr)
    if (errmgr%has_error_occurred()) return

    call eigen(jac, vals, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Test for stability
    rst = .true.
    do i = 1, neqn
        if (real(vals(i)) > 0.0d0) then
            rst = .false.
            return
        end if
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("is_stable_equilibrium", trim(errmsg), &
        DIFFEQ_MEMORY_ALLOCATION_ERROR)
    return

    ! Eigenvalue Array Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The eigenvalue array was expected to be of size ", &
        neqn, ", but was found to be of size ", size(lambda), "."
    call errmgr%report_error("is_stable_equilibrium", trim(errmsg), &
        DIFFEQ_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end function

! ------------------------------------------------------------------------------
end module