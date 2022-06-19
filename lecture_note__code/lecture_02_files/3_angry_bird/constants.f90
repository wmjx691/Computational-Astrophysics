!---------------------------------------------------
! The constants module
!
module constants

    implicit none
    
    real, parameter :: c = 2.99792458e8 ! m/s
    real, parameter :: g = 9.8          ! m/s/s
    real, parameter :: pi = 4.0*atan(1.0)

    contains

        subroutine show_constants()
            implicit none
            print *, "g  =", g
            print *, "pi =", pi

        end subroutine show_constants

end module constants
