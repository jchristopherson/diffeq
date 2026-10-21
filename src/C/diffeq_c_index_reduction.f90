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

module diffeq_c_index_reduction
    use iso_c_binding
    use iso_fortran_env
    use diffeq_c_function_pointers
    use diffeq_index_reduction
    use diffeq_c_errors
    use diffeq_base
    implicit none
    private
    public :: c_ggl_container

    type, bind(C) :: c_ggl_container
        !! C-interoperable description of a holonomically constrained
        !! mechanical model to be solved with the stabilized index-1
        !! Gear-Gupta-Leimkuhler (GGL) formulation of [[ggl_container]].
        !!
        !! For N coordinates and M constraints the integrated state is
        !! \[
        !! z=[q_N,v_N,\lambda_M,\mu_M]^T,
        !! \]
        !! so the state count supplied to the DAE integrators must equal
        !! \(2(N+M)\), and the initial conditions must use this ordering.
        !! The model is
        !! \[
        !! \dot q = v + G^T\mu, \qquad M(q)\dot v = F(x,q,v) - G^T\lambda,
        !! \]
        !! \[
        !! 0 = G\dot q + \alpha g, \qquad
        !! 0 = DG(q)[\dot q]v + GM^{-1}(F-G^T\lambda) + \beta Gv,
        !! \]
        !! where \(g(q)\) are the position constraints and
        !! \(G=\partial g/\partial q\). Every callback must accept perturbed,
        !! off-manifold states, and the constraint callbacks must not depend
        !! on the independent variable.
        type(c_funptr) :: force
            !! Required [[c_ode]] callback computing the applied force
            !! \(F(x,q,v)\), excluding multiplier forces. It is called with
            !! n = 2N and receives the 2N-element vector [q, v], but must
            !! write only the N-element force vector.
        type(c_funptr) :: mass
            !! Optional [[c_ode_mass_matrix]] callback receiving the
            !! N-element coordinate vector q and writing the N-by-N
            !! symmetric positive-definite physical mass matrix. Set to
            !! C_NULL_FUNPTR to use the identity matrix. Set the
            !! mass_matrix_constant setting to false whenever this matrix
            !! varies with x or q.
        type(c_funptr) :: constraints
            !! Required [[c_dae_algebraic_constraints]] callback receiving
            !! the N-element coordinate vector q and writing the M position
            !! constraint residuals \(g(q)\).
        type(c_funptr) :: gradient
            !! Required [[c_ode_jacobian]] callback receiving the N-element
            !! coordinate vector q and writing the M-by-N constraint
            !! Jacobian \(G=\partial g/\partial q\) in column-major order.
            !! The rows of \(G\) must be linearly independent.
        type(c_funptr) :: curvature
            !! Required [[c_dae_constraint_curvature]] callback computing the
            !! M curvature terms \(DG(q)[\dot q]v\). Return zeros for linear
            !! constraints.
        real(c_double) :: alpha
            !! Nonnegative, finite position-constraint decay rate, in inverse
            !! independent-variable units. Zero disables position damping.
        real(c_double) :: beta
            !! Nonnegative, finite velocity-constraint decay rate, in inverse
            !! independent-variable units. Zero disables velocity damping.
        integer(c_int) :: constraint_count
            !! The number M of independent constraints, with
            !! \(1 \le M \le N\).
    end type

end module