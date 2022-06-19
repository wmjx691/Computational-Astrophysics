module Simulation_data
    implicit none

    integer,parameter :: N=9    ! number of stars

    type Star
        integer :: id
        real :: mass
        real :: x
        real :: y
        real :: vx
        real :: vy
        real :: ax
        real :: ay
    end type Star

    type(Star), dimension(N) :: stars
    type(Star), dimension(N) :: stars2
    type(Star), dimension(N) :: stars3
    type(Star), dimension(N) :: stars4

    real :: separation, period

end module Simulation_data
