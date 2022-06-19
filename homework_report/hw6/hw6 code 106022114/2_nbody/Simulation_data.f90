module Simulation_data
    implicit none

    integer,parameter :: N=1    ! number of stars

    type Star
        ! integer :: id
        ! real :: mass
        real*16 :: r
        real*16 :: theta
        real*16 :: r1
        real*16 :: theta1
        real*16 :: r2
        real*16 :: theta2
        real*16 :: l
        real*16 :: energy
        real*16 :: J
        real*16 :: Jin
    end type Star

    type(Star) :: stars
    type(Star) :: stars2
    type(Star) :: stars3
    type(Star) :: stars4

    real :: separation, period

end module Simulation_data
