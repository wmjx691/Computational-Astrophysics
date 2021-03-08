!---------------------------------------------------
! The physics module
!
module physics
    implicit none
    contains
        subroutine update(dt,x0,y0,vx0,vy0,x,y,vx,vy)
            use constants, only : g
            implicit none
            real, intent(in)  :: dt,x0, y0, vx0, vy0
            real, intent(out) :: x,y,vx,vy
            return
        end subroutine update

end module physics

