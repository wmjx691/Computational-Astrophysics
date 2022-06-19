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
    real             :: FL,FR  !flux
    
    call bounary(u)
    uold = u
    do i = istart, iend
      ! finite difference method
      !u(i) = uold(i) - c*dt/dx*(uold(i)-uold(i-1))
      
      !finite volume method
      call flux(i,dt,FL,FR)
      u(i) = uold(i) - dt/dx*(FR-FL)
    enddo
    
end subroutine update


subroutine flux(i, dt,FL,FR)
    use Simulation_data
    implicit none
    integer,intent(in) :: i
    real, intent(in)   :: dt
    real, intent(out)  :: FL
    real, intent(out)  :: FR
    
    real  ::  sig
    
    ! arithmetic average (unstable)
    !FL = 0.5*c*(uold(i-1)+uold(i))
    !FR = 0.5*c*(uold(i+1)+uold(i))
    
    
    !The Lax-Fried
    !FL = 0.5*c*(uold(i-1)+uold(i)) - 0.5*dx/dt*(uold(i)-uold(i-1))
    !FR = 0.5*c*(uold(i+1)+uold(i)) - 0.5*dx/dt*(uold(i+1)-uold(i))
    
    
    !The upwind method
    !FL = uold(i-1)*c
    !FR = uold(i)*c
    
    
    !The Lax-Wendroff Method
    !FL = 0.5*c*(uold(i-1)+uold(i)) - 0.5*dx/dt*c**2*(uold(i)-uold(i-1))
    !FR = 0.5*c*(uold(i+1)+uold(i)) - 0.5*dx/dt*c**2*(uold(i+1)-uold(i))
    
    call reconstruct(dx,uold(i-2),uold(i-1),uold(i),sig)
    FL = c*uold(i-1)+0.5*c*(dx - c*dt)*sig    
    call reconstruct(dx,uold(i-1),uold(i),uold(i+1),sig)
    FR = c*uold(i)  +0.5*c*(dx - c*dt)*sig  
    return

end subroutine flux


subroutine reconstruct(dx,l,m,r,sig)
    implicit none
    real, intent(in)   :: dx,l,m,r
    real, intent(out)  :: sig
    
    real  :: a,b
    
    a = (m-l)/dx
    b = (r-m)/dx
    
    !mimod limiter
    if(a*b .gt. 0.0) then
        if (abs(a) .lt. abs(b)) then
            sig = a
        else
            sig = b
        endif
    else
        sig = 0.0
    endif

    return
end subroutine reconstruct