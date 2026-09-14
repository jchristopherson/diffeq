!> @brief This module contains several ODE solvers and associated types.
module diffeq
    use iso_fortran_env
    use diffeq_base
    use diffeq_runge_kutta
    use diffeq_multistep
    use diffeq_bdf
    use diffeq_pece
    use diffeq_rosenbrock
    use diffeq_kennedy_carpenter
end module