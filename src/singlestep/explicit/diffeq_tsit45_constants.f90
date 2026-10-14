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

module diffeq_tsit45_constants
    !! Coefficients for the Tsitouras 5/4 embedded Runge--Kutta method.
    use iso_fortran_env, only : real64
    implicit none
    private
    public :: tsit_c, tsit_a, tsit_b, tsit_e

    real(real64), parameter :: tsit_c(7) = [ &
        0.0d0, 0.161d0, 0.327d0, 0.9d0, 0.980025540904509685729810680793d0, &
        1.0d0, 1.0d0 &
    ]

    real(real64), parameter :: tsit_a(7,7) = reshape([ &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.161d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        -0.008480655492356989d0, 0.335480655492356989d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        2.8971530571054935d0, -6.359448489975075d0, 4.3622954328695815d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        5.325864828439257d0, -11.748883564062828d0, 7.495539342889836d0, -0.09249506636175525d0, 0.0d0, 0.0d0, 0.0d0, &
        5.861455442946420d0, -12.92096931784711d0, 8.159367898576159d0, -0.071584973281401d0, &
        -0.02826905039406838d0, 0.0d0, 0.0d0, &
        0.09646076681806523d0, 0.01d0, 0.4798896504144996d0, 1.379008574103742d0, -3.290069515436081d0, 2.324710524099774d0, 0.0d0 &
    ], [7,7], order = [2,1])

    real(real64), parameter :: tsit_b(7) = [ &
        0.09646076681806523d0, 0.01d0, 0.4798896504144996d0, 1.379008574103742d0, &
        -3.290069515436081d0, 2.324710524099774d0, 0.0d0 &
    ]

    real(real64), parameter :: tsit_e(7) = [ &
        -0.001780011052226d0, -0.0008164344596567469d0, 0.007880878010261995d0, &
        -0.1447110071732629d0, 0.5823571654525552d0, -0.45808210592918697d0, &
        0.0151515151515151515151515151515d0 &
    ]
end module diffeq_tsit45_constants
