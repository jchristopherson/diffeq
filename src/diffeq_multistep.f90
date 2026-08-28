module diffeq_multistep
    !! Base definitions for modern multistep ODE solvers.
    !!
    !! Concrete multistep methods will be added here without relying on
    !! legacy fixed-form solver implementations.
    use diffeq_base, only : ode_integrator
    implicit none
    private
    public :: multistep_integrator

    type, abstract, extends(ode_integrator) :: multistep_integrator
        !! Abstract base type for variable-step multistep integrators.
    end type

end module diffeq_multistep
