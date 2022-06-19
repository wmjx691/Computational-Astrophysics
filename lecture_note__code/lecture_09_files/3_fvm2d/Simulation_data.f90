module Simulation_data
        implicit none
        integer, parameter  ::  imax   = 120  !number of points in x-direction
        integer, parameter  ::  ibuf   = 1    !number of ghost zones for B.C.
        integer, parameter  ::  istart = 1
        integer, parameter  ::  iend   = imax

        integer, parameter  ::  jmax   = 120  !number of points in y-direction
        integer, parameter  ::  jbuf   = 1    !number of ghost zones for B.C.
        integer, parameter  ::  jstart = 1
        integer, parameter  ::  jend   = jmax
        
        real, parameter     ::  cx   = 1.0
        real, parameter     ::  cy   = 1.0
        real, parameter     ::  xmin = 0.0
        real, parameter     ::  xmax = 1.0
        real, parameter     ::  ymin = 0.0
        real, parameter     ::  ymax = 1.0
        real, parameter     ::  tend = 1.0

        real, parameter     ::  cfl  = 0.4
        real, save          ::  dx, dy

        real, dimension(istart-ibuf:iend+ibuf,jstart-jbuf:jend+jbuf), save  :: u, uold
        real, dimension(istart-ibuf:iend+ibuf), save  :: x
        real, dimension(jstart-jbuf:jend+jbuf), save  :: y

        integer, parameter  ::  io_interval = 5

end module Simulation_data
