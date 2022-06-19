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

            ! x0,y0,vx0,vy0 are input positions and velocities
            ! x,y,vx,vy are values at t+dt step

            x = x0 + vx0*dt
            y = y0 + vy0*dt
            vx = vx0
            vy = vy0 - g*dt
            return
        end subroutine update

end module physics

