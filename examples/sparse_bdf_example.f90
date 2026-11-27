! This file is part of diffeq.
!
! diffeq is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! diffeq is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with diffeq. If not, see <https://www.gnu.org/licenses/>.

program sparse_bdf_example
    !! Solve the one-dimensional heat equation using a sparse BDF Jacobian.
    !!
    !! Centered differences on 100 interior points produce a stiff system
    !! with a tridiagonal Jacobian. Both boundaries are held at zero. The
    !! initial condition combines two discrete eigenmodes; their known decay
    !! rates provide an independent check of the solution at each output time.
    use iso_fortran_env, only : real64
    use diffeq
    use fplot_core
    use linalg, only : csr_matrix, create_csr_matrix
    implicit none

    integer, parameter :: n = 1000, noutput = 5, high_mode = 15
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: dx = 1.0d0 / (n + 1)
    real(real64), parameter :: diffusivity = 1.0d0
    real(real64), parameter :: rate = diffusivity / dx**2
    real(real64), parameter :: times(noutput) = [0.0d0, 0.002d0, 0.01d0, 0.05d0, 0.1d0]
    real(real64) :: initial(n), position(n+2), profile(n+2), expected(n), max_error, lambda1, lambda15
    real(real64), allocatable :: solution(:,:)
    type(bdf) :: integrator
    type(ode_container) :: model
    type(plot_2d) :: plt
    class(terminal), pointer :: term
    type(plot_data_2d) :: curves(noutput)
    class(legend), pointer :: key
    character(len=32) :: label
    integer :: i, j, jacobian_calls

    ! The sparse Jacobian callback selects BDF's sparse Newton path.
    model%fcn => heat_rhs
    model%sparse_jacobian => heat_jacobian
    jacobian_calls = 0
    do i = 1, n+2
        position(i) = real(i-1, real64) * dx
    end do
    initial = sin(pi * position(2:n+1)) + 0.5d0 * sin(high_mode * pi * position(2:n+1))

    ! Integrate
    call integrator%solve(model, times, initial)
    solution = integrator%get_solution()

    ! Exact eigenvalues of the discrete Laplacian, not the continuous PDE.
    lambda1 = -4.0d0 * rate * sin(pi / (2.0d0 * (n+1)))**2
    lambda15 = -4.0d0 * rate * sin(high_mode * pi / (2.0d0 * (n+1)))**2
    max_error = 0.0d0
    do j = 1, noutput
        expected = exp(lambda1 * times(j)) * sin(pi * position(2:n+1)) + &
            0.5d0 * exp(lambda15 * times(j)) * sin(high_mode * pi * position(2:n+1))
        max_error = max(max_error, maxval(abs(solution(j,2:n+1) - expected)))
    end do
    print "(A,I0,A,I0)", "Sparse system: ", n, " equations, Jacobian nonzeros: ", 3*n-2
    print "(A,I0)", "Sparse Jacobian evaluations: ", jacobian_calls
    print "(A,ES12.4)", "Maximum error against discrete exact solution: ", max_error
    if (jacobian_calls == 0) error stop "Sparse Jacobian callback was not invoked"
    if (max_error > 1.0d-4) error stop "Sparse BDF solution failed the accuracy check"

    ! Plotting
    call plt%initialize()
    term => plt%get_terminal()
    call term%set_window_width(800)
    call plt%set_x_axis_title("Position")
    call plt%set_y_axis_title("Temperature")
    key => plt%get_legend()
    call key%set_is_visible(.true.)
    call key%set_draw_border(.false.)
    call key%set_draw_inside_axes(.false.)
    call key%set_vertical_position(LEGEND_CENTER)
    do j = 1, noutput
        profile(1) = 0.0d0
        profile(2:n+1) = solution(j,2:n+1)
        profile(n+2) = 0.0d0
        call curves(j)%define_data(position, profile)
        write(label, '("t = ",F5.3)') times(j)
        call curves(j)%set_name(trim(label))
        call curves(j)%set_line_width(2.0)
        call plt%push(curves(j))
    end do
    call plt%draw()

contains
    subroutine heat_rhs(t, y, dydt, args)
        !! Apply the finite-difference Laplacian, including zero boundaries.
        real(real64), intent(in) :: t, y(:)
        real(real64), intent(out) :: dydt(:)
        class(*), intent(in), optional :: args
        integer :: k

        dydt(1) = rate * (-2.0d0 * y(1) + y(2))
        do k = 2, n-1
            dydt(k) = rate * (y(k-1) - 2.0d0 * y(k) + y(k+1))
        end do
        dydt(n) = rate * (y(n-1) - 2.0d0 * y(n))
    end subroutine

    subroutine heat_jacobian(t, y, jac, args)
        !! Supply the tridiagonal Jacobian directly as a CSR matrix.
        real(real64), intent(in) :: t, y(:)
        type(csr_matrix), intent(out) :: jac
        class(*), intent(in), optional :: args
        integer :: rows(3*n-2), columns(3*n-2), k, entry
        real(real64) :: values(3*n-2)

        entry = 0
        do k = 1, n
            if (k > 1) then
                entry = entry + 1
                rows(entry) = k
                columns(entry) = k-1
                values(entry) = rate
            end if
            entry = entry + 1
            rows(entry) = k
            columns(entry) = k
            values(entry) = -2.0d0 * rate
            if (k < n) then
                entry = entry + 1
                rows(entry) = k
                columns(entry) = k+1
                values(entry) = rate
            end if
        end do
        jacobian_calls = jacobian_calls + 1
        jac = create_csr_matrix(n, n, rows, columns, values)
    end subroutine
end program sparse_bdf_example