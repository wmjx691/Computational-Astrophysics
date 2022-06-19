!---------------------------------------------------
! The main program
!
program nbody

    use constants
    use Simulation_data
    use physics, only : update, initial_from_model
    use IO, only : record

    implicit none

    integer :: interval, noutput, step
    real, parameter :: elipson=0.02
    real, parameter :: omegap=1.5
    real, parameter :: vc=1.0,vr=0.0
    integer, parameter :: m=2
    real :: dt, time, tmax

    step     = 0
    dt       = 0.001    ! sec
    time     = 0.0           ! sec
    tmax     = 200    ! sec
!    tmax     = 1.414213*8.0*3.1415926    ! sec 
    noutput  = 1000
    interval = (tmax/dt)/noutput

    call initial_from_model(elipson,omegap,vc,vr,m)

    do while (time .le. tmax)
        !call update(time,dt)

        call update(time,elipson,omegap,vc,m,dt)
        time = time + dt

        if (mod(step,interval) .eq. 0) then
            call record(time)
        endif
        step = step + 1
        print *, time, stars%r, stars%theta
    end do

end program nbody

