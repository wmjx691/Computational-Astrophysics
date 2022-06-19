module Simulation_data
    implicit none

    integer, parameter :: nx=2000
!
    real, dimension(nx) :: fm
    real, dimension(nx) :: dm
    real, dimension(nx) :: dm12
    real, dimension(nx) :: r
    real, dimension(nx) :: dr
    real, dimension(nx) :: dr12
    real, dimension(nx) :: rold
    real, dimension(nx) :: v
    real, dimension(nx) :: a
    real, dimension(nx) :: at12
    real, dimension(nx) :: ak12
    real, dimension(nx) :: u
    real, dimension(nx) :: p
    real, dimension(nx) :: rho
    real, dimension(nx) :: eps
    real, dimension(nx) :: w
    real, dimension(nx) :: wt12
    real, dimension(nx) :: aux

    logical, save ::  cartesian

    real, save :: pi,gamma,t,dt,dt12,tmax,q,fmratio,&
                  cflfactor,dfactor,efac
    integer, save :: il, istep

end module Simulation_data
