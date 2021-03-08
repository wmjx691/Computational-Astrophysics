program angry_bird

    implicit none

    real :: angle, velocity

    real :: dt, time, velx, vely, posx, posy
    real, parameter :: g = 9.8  ! m/s/s
    real :: pi = 4.0*atan(1.0)

    real :: analy_x, analy_y

    ! Input parameters
    angle    = 60.0                ! degree
    angle    =                     ! change to rad
    velocity = 30.0 ! m/s 
    velx = 
    vely = 

    print *, "Choose a time step in second"
    read *, dt ! sec

    open(unit=11, file="output.txt")


    ! initial conditions
    time = 
    posx = 
    posy = 
    analy_x =   ! analytical solution
    analy_y = 

    write(11,100) "#", "time", "posx", "posy","anax"," anay"
    write(11,200) time, posx, posy, analy_x, analy_y


    do while (posy .ge. 0.0)

        !
        ! In this example we use a first order scheme (Euler method)
        ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
        ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
        !
        ! the same approximation can be applied to dv/dt = a
        !


        ! update vel and pos for one step


        ! update time


        ! compare it with the analytical solution
        analy_x = (velocity * cos(angle)) * time
        analy_y = velocity * sin(angle) * time - 0.5*g*time**2

        print *, time, posx, posy, analy_x, analy_y
        write(11,200) time, posx, posy, analy_x, analy_y
    end do

    close(11)
100 format(a2, 5a24)
200 format(2x, 5e24.14)

end program

