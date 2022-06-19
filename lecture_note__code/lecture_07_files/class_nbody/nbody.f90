!---------------------------------------------------
! The main program
!
program nbody

    use constants
    use Simulation_data
    use physics, only : initial_from_model, update
    use IO, only : record

    implicit none

    integer :: i, interval, noutput, step
    real :: dt, time, tmax

    step     = 0
    dt       = 1.e-3 * yr    ! sec
    time     = 0.0           ! sec
    tmax     = 100.0 * yr    ! sec 
    noutput  = 1000
    interval = (tmax/dt)/noutput

    call initial_from_model()

    do while (time .le. tmax)
        !call update(time,dt)

        call update(time,dt)
        time = time + dt

        if (mod(step,interval) .eq. 0) then
            call record(time)
        endif
        step = step + 1
        !print *, time, stars(1)%x, stars(1)%y, stars(2)%x, stars(2)%y
    end do

end program nbody

