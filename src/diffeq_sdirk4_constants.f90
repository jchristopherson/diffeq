module diffeq_sdirk4_constants
    use iso_fortran_env
    implicit none
    
    ! Model Parameters (ESDIRK4(3)6L[2]SA, Table 16, pg 90)
    ! https://ntrs.nasa.gov/api/citations/20160005923/downloads/20160005923.pdf
    real(real64), parameter :: gamma = 0.25d0
    real(real64), parameter :: c2 = 0.5d0
    real(real64), parameter :: c3 = (2.0d0 - sqrt(2.0d0)) / 4.0d0
    real(real64), parameter :: c4 = 5.0d0 / 8.0d0
    real(real64), parameter :: c5 = 26.0d0 / 25.0d0
    real(real64), parameter :: c6 = 1.0d0
    real(real64), parameter :: a21 = gamma
    real(real64), parameter :: a22 = gamma
    real(real64), parameter :: a32 = (1.0d0 - sqrt(2.0d0)) / 8.0d0
    real(real64), parameter :: a31 = (c3 - a32 - gamma)
    real(real64), parameter :: a33 = gamma
    real(real64), parameter :: a43 = 7.0d0 * (1.0d0 + sqrt(2.0d0)) / 32.0d0
    real(real64), parameter :: a42 = (5.0d0 - 7.0d0 * sqrt(2.0d0)) / 64.0d0
    real(real64), parameter :: a41 = c4 - a42 - a43 - gamma
    real(real64), parameter :: a44 = gamma
    real(real64), parameter :: a54 = &
        166.0d0 * (-97.0d0 + 376.0d0 * sqrt(2.0d0)) / 1.09375d5
    real(real64), parameter :: a53 = &
        (5.06605d5 + 1.32109d5 * sqrt(2.0d0)) / 4.375d5
    real(real64), parameter :: a52 = &
        (-1.3796d4 - 5.4539d4 * sqrt(2.0d0)) / 1.25d5
    real(real64), parameter :: a51 = c5 - a52 - a53 - a54 - gamma
    real(real64), parameter :: a55 = gamma
    real(real64), parameter :: b5 = &
        -15.625d3 * (97.0d0 + 376.0d0 * sqrt(2.0d0)) / 9.0749876d7
    real(real64), parameter :: b4 = &
        -16.0d0 * (-2.2922d4 + 3.525d3 * sqrt(2.0d0)) / 5.71953d5
    real(real64), parameter :: b3 = &
        47.0d0 * (-267.0d0 + 1.783d3 * sqrt(2.0d0)) / 2.73343d5
    real(real64), parameter :: b2 = (1.181d3 - 9.87d2 * sqrt(2.0d0)) / 1.3782d4
    real(real64), parameter :: b1 = 1.0d0 - b2 - b3 - b4 - b5 - gamma
    real(real64), parameter :: b6 = gamma
    real(real64), parameter :: b2a = -4.80923228411d11 / 4.982971448372d12
    real(real64), parameter :: b3a = 6.709447293961d12 / 1.2833189095359d13
    real(real64), parameter :: b4a = 3.513175791894d12 / 6.748737351361d12
    real(real64), parameter :: b5a = -4.98863281070d11 / 6.042575550617d12
    real(real64), parameter :: b6a = 2.077005547802d12 / 8.945017530137d12
    real(real64), parameter :: b1a = 1.0d0 - b2a - b3a - b4a - b5a - b6a

    ! Dense Output Coefficients
    real(real64), parameter :: bs11 = 11963910384665.0d0 / 12483345430363.0d0
    real(real64), parameter :: bs12 = 11963910384665.0d0 / 12483345430363.0d0
    real(real64), parameter :: bs13 = -28603264624.0d0 / 1970169629981.0d0
    real(real64), parameter :: bs14 = -3524425447183.0d0 / 2683177070205.0d0
    real(real64), parameter :: bs15 = -17173522440186.0d0 / 10195024317061.0d0
    real(real64), parameter :: bs16 = 27308879169709.0d0 / 13030500014233.0d0
    real(real64), parameter :: bs21 = -69996760330788.0d0 / 18526599551455.0d0
    real(real64), parameter :: bs22 = -69996760330788.0d0 / 18526599551455.0d0
    real(real64), parameter :: bs23 = 102610171905103.0d0 / 26266659717953.0d0
    real(real64), parameter :: bs24 = 74957623907620.0d0 / 12279805097313.0d0
    real(real64), parameter :: bs25 = 113853199235633.0d0 / 9983266320290.0d0
    real(real64), parameter :: bs26 = -84229392543950.0d0 / 6077740599399.0d0
    real(real64), parameter :: bs31 = 32473635429419.0d0 / 7030701510665.0d0
    real(real64), parameter :: bs32 = 32473635429419.0d0 / 7030701510665.0d0
    real(real64), parameter :: bs33 = -38866317253841.0d0 / 6249835826165.0d0
    real(real64), parameter :: bs34 = -26705717223886.0d0 / 4265677133337.0d0
    real(real64), parameter :: bs35 = -121105382143155.0d0 / 6658412667527.0d0
    real(real64), parameter :: bs36 = 1102028547503824.0d0 / 51424476870755.0d0
    real(real64), parameter :: bs41 = -14668528638623.0d0 / 8083464301755.0d0
    real(real64), parameter :: bs42 = -14668528638623.0d0 / 8083464301755.0d0
    real(real64), parameter :: bs43 = 21103455885091.0d0 / 7774428730952.0d0
    real(real64), parameter :: bs44 = 30155591475533.0d0 / 15293695940061.0d0
    real(real64), parameter :: bs45 = 119853375102088.0d0 / 14336240079991.0d0
    real(real64), parameter :: bs46 = -63602213973224.0d0 / 6753880425717.0d0
end module