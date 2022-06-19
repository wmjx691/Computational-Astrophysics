module Simulation_data
    implicit none

    integer,parameter :: N=11    ! number of stars

    type Star
        integer :: id  
        real :: mass
        real :: x
        real :: y
        real :: vx
        real :: vy
        real :: ax
        real :: ay
        real :: xt
        real :: yt
        real :: vxt
        real :: vyt
        real :: axt
        real :: ayt

    end type Star

    type(Star), dimension(N) :: stars
    real :: separation
    real :: period

end module Simulation_data
