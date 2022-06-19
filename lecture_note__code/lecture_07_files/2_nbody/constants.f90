!---------------------------------------------------
! The constants module
!
module constants

    implicit none

    real, parameter :: G= 6.67428e-8    !
    real, parameter :: pi = 4.0*atan(1.0)
    real, parameter :: msun = 1.989e33
    real, parameter :: au = 1.49598e13
    real, parameter :: yr = 365.*86400.0

    contains

        subroutine show_constants()
            implicit none
            print *, "G  =", G
            print *, "pi =", pi

        end subroutine show_constants

end module constants

