module Simulation_data
        implicit none
        integer, parameter  ::  imax   = 20   !number of points in x-direction
        integer, parameter  ::  ibuf   = 1    !number of ghost zones for B.C.
        integer, parameter  ::  istart = 0
        integer, parameter  ::  iend   = imax

        real, parameter     ::  c    = 1.0
        real, parameter     ::  xmin = 0.0
        real, parameter     ::  xmax = 1.0
        real, parameter     ::  tend = 0.06

        real, parameter     ::  cfl  = 0.52
        real, save          ::  dx

        real, dimension(istart-ibuf:iend+ibuf), save  :: u, uold, x

        integer, parameter  ::  io_interval = 10

end module Simulation_data
