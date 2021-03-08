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

    do while (posy .ge. 0.0)
        call update(dt,posx,posy,velx,vely,posx,posy,velx,vely)
        time = time + dt
        print *, time, posx, posy, velx, vely
    end do

end program angry_bird

