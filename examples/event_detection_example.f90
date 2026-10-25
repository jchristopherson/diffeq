! This file is part of diffeq.
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

program example
    use iso_fortran_env
    use diffeq
    use diffeq_models
    use fplot_core
    implicit none

    ! Model Parameters
    real(real64), parameter :: m = 1.0d0
    real(real64), parameter :: g = 9.81d0
    real(real64), parameter :: k = 0.7d0
    real(real64), parameter :: initial_height = 5.0d0
    real(real64), parameter :: initial_velocity = 0.0d0

    ! Local Variables
    type(ode_container) :: model
    type(tsitouras_54) :: integrator
    integer :: n
    real(real64) :: ic(2), t(2), v
    real(real64), allocatable, dimension(:) :: kinetic_energy, &
        potential_energy, total_energy
    real(real64), allocatable, dimension(:,:) :: sol

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    class(legend), pointer :: lgnd
    class(plot_axis), pointer :: xAxis
    
    ! Set up the model and the event detection
    model%fcn => ball
    model%event_function => ball_ground_constraint
    call model%set_event_equation_count(1)  ! we only have 1 event equation

    ! Define the initial conditions and time vector
    t = [0.0d0, 1.0d1]
    ic = [initial_height, initial_velocity]

    ! Compute the solution.  The solver will return at the location of an event
    call integrator%set_maximum_step_size(1.0d-2) ! not necessary, but here to illustrate a smooth solution plot
    call integrator%solve(model, t, ic)

    ! Now that the ball has hit the ground, we need to reset the initial 
    ! conditions to continue the bouncing.  We'll need the final velocity of
    ! the ball as it impacted to properly reset the initial conditions.  Also
    ! notice that we do not have to store the initial conditions for this solve
    ! in the solution buffer as the impact point from the previous solve is
    ! already stored.  This can be accomplished via the optional store_iv
    ! option in the solver.
    sol = integrator%get_solution()
    n = size(sol, 1)
    v = sol(n,3)
    ic = [0.0d0, -k * v]    ! velocity is in the opposite direction; hence, the sign change
    t(1) = sol(n, 1)        ! also reset the integration start time
    call integrator%solve(model, t, ic, store_iv = .false.)
    
    ! Let's repeat that for another bounce
    sol = integrator%get_solution()
    n = size(sol, 1)
    v = sol(n, 3)
    ic = [0.0d0, -k * v]
    t(1) = sol(n, 1)
    call integrator%solve(model, t, ic, store_iv = .false.)
    sol = integrator%get_solution()

    ! For illustration purposes, let's compute the kinetic energy of the ball
    kinetic_energy = 0.5d0 * m * sol(:,3)**2
    potential_energy = m * g * sol(:,2)
    total_energy = kinetic_energy + potential_energy

    ! Plot the solution
    call plt%initialize(1, 2, width = 1200)

    call plt1%initialize()
    xAxis => plt1%get_x_axis()
    call xAxis%set_zero_axis(.true.)
    call plt1%set_x_axis_title("t")
    call plt1%set_y_axis_title("y(t)")
    call plt1%push(sol(:,1), sol(:,2))
    call plt%set(1, 1, plt1)

    call plt2%initialize()
    xAxis => plt2%get_x_axis()
    lgnd => plt2%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call xAxis%set_zero_axis(.true.)
    call plt2%set_x_axis_title("t")
    call plt2%set_y_axis_title("E")
    call plt2%push(sol(:,1), kinetic_energy, name = "Kinetic Energy")
    call plt2%push(sol(:,1), potential_energy, ls = LINE_DASHED, name = "Potential Energy")
    call plt2%push(sol(:,1), total_energy, lw = 1.2, ls = LINE_DASH_DOTTED, &
        lc = CLR_BLACK, name = "Total Energy")
    call plt%set(1, 2, plt2)

    call plt%draw()
end program