!===================================================================
!
! The simpliest first order finite-difference method for calculating 
!
! the advection equation: u_t = - c u_x
!
! 2019/04/23  Instittue of Astronomy, NTHU
! 2020/05/13  Instittue of Astronomy, NTHU
! 
!===================================================================
program advection

    implicit none

    call initial()

    call evolution()

    call finalize()

end program advection
