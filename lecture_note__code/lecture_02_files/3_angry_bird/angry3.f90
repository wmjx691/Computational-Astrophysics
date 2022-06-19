!---------------------------------------------------
! The main program
!
program angry_bird

    use constants, only : show_constants, g ,pi, c
    use physics, only : update

    implicit none

    real :: angle, velocity

    real :: dt, time, velx, vely, posx, posy

    !TODO
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

