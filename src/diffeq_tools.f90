!
module diffeq_tools
    use iso_fortran_env
    use diffeq
    use ferror
    use nonlin_core
    use nonlin_solve
    implicit none
    private
    public :: find_equilibrium_points

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
    call solverptr%solve(container, rst(:,1), rst(:,2))

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
! TO DO: Test for stability of equilibrium points.
!
! Assume: dy/dx = f(x), and x* is the equilibrium point
! f'(x*) < 0 then x* is a stable equilibrium point
! f'(x*) > 0 then x* is an unstable equilibrium point
!
! For systems of equations, compute the Jacobian matrix and consider the
! eigenvalues.  If all eigenvalues are negative, the equilibrium point is
! stable.  If all eigenvalues are positive, the equilibrium point is unstable.
! If the eigenvalues are mixed, it is an unstable saddle point

! ------------------------------------------------------------------------------
end module