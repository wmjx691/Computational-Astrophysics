subroutine evolution
    use IO, only: output
    use Simulation_data
    implicit none

    integer  ::  n
    real     ::  dt, time

    n   = 0
    time = 0.0
    
    dt = min(abs(dx/cx),abs(dy/cy))*cfl

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
    integer          :: i,j
    real             :: FL,FR  !flux
    
    !do x-direction
    call boundary_x(u)
    uold = u
    do j = jstart, jend
        do i = istart, iend      
          !finite volume method
          call flux_x(i,j,dt,FL,FR)
          u(i,j) = uold(i,j) - dt/dx*(FR-FL)
        enddo
    enddo
    
    !do y-direction
    call boundary_y(u)
    uold = u
    do j = jstart, jend
        do i = istart, iend      
          !finite volume method
          call flux_y(i,j,dt,FL,FR)
          u(i,j) = uold(i,j) - dt/dy*(FR-FL)
        enddo
    enddo
    
end subroutine update


subroutine flux_x(i,j,dt,FL,FR)
    use Simulation_data
    implicit none
    integer,intent(in) :: i,j
    real, intent(in)   :: dt
    real, intent(out)  :: FL
    real, intent(out)  :: FR
    
    real  ::  sig
    
    call reconstruct(dx,uold(i-2,j),uold(i-1,j),uold(i,j),sig)
    FL = cx*uold(i-1,j)+0.5*cx*(dx - cx*dt)*sig    
    call reconstruct(dx,uold(i-1,j),uold(i,j),uold(i+1,j),sig)
    FR = cx*uold(i,j)  +0.5*cx*(dx - cx*dt)*sig  
    return

end subroutine flux_x


subroutine flux_y(i,j,dt,FL,FR)
    use Simulation_data
    implicit none
    integer,intent(in) :: i,j
    real, intent(in)   :: dt
    real, intent(out)  :: FL
    real, intent(out)  :: FR
    
    real  ::  sig
    
    call reconstruct(dy,uold(i,j-2),uold(i,j-1),uold(i,j),sig)
    FL = cy*uold(i,j-1)+0.5*cy*(dy - cy*dt)*sig    
    call reconstruct(dy,uold(i,j-1),uold(i,j),uold(i,j+1),sig)
    FR = cy*uold(i,j)  +0.5*cy*(dy - cy*dt)*sig  
    return

end subroutine flux_y


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