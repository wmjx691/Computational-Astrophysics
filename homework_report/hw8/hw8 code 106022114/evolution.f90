subroutine evolution
    use IO, only: output
    use Simulation_data
    implicit none

    integer  ::  n
    real     ::  dt, time

    n   = 0
    time = 0.0
    
    dt = 0.0013

    ! main loop

    do while(time .le. tend)

        if (mod(n,io_interval) .eq. 0) then   !record it every io_intercal
            print *, "n=",n, "Time = ", time
            call output(n,time)
        endif    
        
        ! update u
        call update(time, dt)

        n = n + 1
        time = time + dt

        if (n .eq. 46) then   !record it every io_intercal

            print *, "n=",n, "Time = ",time
            call output(n,time)
        endif

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
      u(i) = uold(i) + (dt/(dx**2))*(uold(i+1)-2.0*uold(i)+uold(i-1))
    enddo

    call bounary(u)
    
end subroutine update