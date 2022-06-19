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

!---------------------------------------------------
! The main program
!
program angry_bird

    use constants, only : show_constants, g ,pi, c
    use physics, only : update

    implicit none

    real :: angle, velocity

    real :: dt, time, velx, vely, posx, posy

    call show_constants()

    ! TODO
    dt = 0.001
    time = 0.0
    posx = 0.0
    posy = 0.0
    velx = 15.0
    vely = 15.0*cos(pi/3) 

    open(unit=11, file="output1.txt")

    do while (posy .ge. 0.0)
        call update(dt,posx,posy,velx,vely,posx,posy,velx,vely)
        time = time + dt
        print *, time, posx, posy, velx, vely
        write(11,*) time, posx, posy, velx, vely
    end do

    close(11)

end program angry_bird

