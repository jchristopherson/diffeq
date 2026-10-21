module diffeq_c_function_pointers
    use iso_c_binding
    implicit none

    abstract interface
        subroutine c_ode(n, x, y, dydx) bind(C)
            !! ODE callback routine.
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
                !! The number of state variables.
            real(c_double), intent(in), value :: x
                !! The current value of the independent variable.
            real(c_double), intent(in) :: y(n)
                !! The current state vector.
            real(c_double), intent(out) :: dydx(n)
                !! The state vector derivatives.
        end subroutine

        subroutine c_ode_jacobian(m, n, x, y, jac, ldj) bind(C)
            !! ODE Jacobian callback routine.
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: m
                !! The number of equations.
            integer(c_int), intent(in), value :: n
                !! The state vector length.
            real(c_double), intent(in), value :: x
                !! The current value of the independent variable.
            real(c_double), intent(in) :: y(n)
                !! The current state vector.
            real(c_double), intent(out) :: jac(ldj,n)
                !! The M-by-N Jacobian matrix.
            integer(c_int), intent(in), value :: ldj
                !! The leading dimension of the Jacobian matrix.
        end subroutine

        subroutine c_ode_mass_matrix(n, x, y, m, ldm) bind(C)
            !! ODE mass matrix callback routine.
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: n
                !! The number of state variables.
            real(c_double), intent(in), value :: x
                !! The current value of the independent variable.
            real(c_double), intent(in) :: y(n)
                !! The current state vector.
            real(c_double), intent(out) :: m(ldm,n)
                !! The mass matrix.
            integer(c_int), intent(in), value :: ldm
                !! The leading dimension of the mass matrix.
        end subroutine
        
        subroutine c_dae_algebraic_constraints(nconstraints, nstates, x, y, g)
            !! Evaluates algebraic constraint residuals \(g(x,y)\).
            !!
            !! For \(N\) input state components and \(M\) constraints,
            !! the callback represents a mapping
            !! \[
            !! g : \mathbb{R}\times\mathbb{R}^{N}
            !! \longrightarrow \mathbb{R}^{M}.
            !! \]
            !! A consistent state satisfies \(g(x,y)=0\). Return the
            !! residuals themselves, not their derivatives, a Boolean result,
            !! or a projected state. A nonzero residual is a valid evaluation
            !! at an inconsistent state, rather than a callback failure.
            !!
            !! The caller supplies the input and output array sizes. The
            !! implementation must interpret the same state ordering as its
            !! model and write every output element. This interface performs
            !! no dimension, rank, or initial-consistency checks. Residual
            !! ordering and physical units are defined by the model; any
            !! subsequent consistency tolerance must respect that scaling.
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: nconstraints
                !! The number of constraint equations.
            integer(c_int), intent(in), value :: nstates
                !! The number of state variables.
            real(c_double), intent(in), value :: x
                !! The current value of the independent variable.
            real(c_double), intent(in) :: y(nstates)
                !! The current state vector.
            real(c_double), intent(out) :: g(nconstraints)
                !! The constraint vector.
        end subroutine

        subroutine c_dae_constraint_curvature(nconstraints, ncoord, x, q, &
            direction, v, curvature)
            !! Evaluates \(DG(q)[\mathrm{direction}]v\), without including
            !! \(G\dot v\). Each component is
            !! \(\mathrm{direction}^T\nabla^2 g_i(q)v\).
            !! Assign every output element. For linear constraints return
            !! zeros. The callback must support off-manifold coordinates and
            !! arbitrary directions used during nonlinear/Jacobian evaluation.
            use iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in), value :: nconstraints
                !! The number of constraint equations.
            integer(c_int), intent(in), value :: ncoord
                !! The number of coordinate variables.
            real(c_double), intent(in), value :: x
                !! The current value of the independent variable.
            real(c_double), intent(in) :: q(ncoord)
                !! The current coordinate vector.
            real(c_double), intent(in) :: direction(ncoord)
                !! The coordinate derivative \(\dot q\), including the 
                !! kinematic multiplier correction.
            real(c_double), intent(in) :: v(ncoord)
                !! The physical velocity components.
            real(c_double), intent(out) :: curvature(nconstraints)
                !! The curvature terms in the order of teh constraint callback.
        end subroutine
    end interface

end module