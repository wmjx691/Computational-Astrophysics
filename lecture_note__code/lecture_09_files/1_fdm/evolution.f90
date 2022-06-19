subroutine evolution
    use IO, only: output
    use Simulation_data
    implicit none

    integer  ::  n
    real     ::  dt, time

    n   = 0
    time = 0.0
    
    dt = abs(dx/c)*cfl

    ! main loop

    do while(time .le. tend)

    if (mod(n,io_interval) .eq. 0) then
        print *, "n=",n, "Time = ", time
        call output(n,time)
        endif    
        ! update u
        call update(time, dt)

        n = n + 1
        time = time + dt
    enddo

end subroutine evolution



subroutine update(time, dt)
    use Simulation_data
    implicit none
    real, intent(in) :: time,dt
    integer          :: i
    
    call bounary(u)
    uold = u
    do i = istart, iend
      u(i) = uold(i) - c*dt/dx*(uold(i)-uold(i-1))
    enddo
    
end subroutine update