!---------------------------------------------------
! The main program
!
program polytrope
    implicit none
!
!    real :: dt, time, vel, pos, m
!
!    dt   = 0.001  
!    time = 0.0
!
!    pos  = 1.0
!    vel  = 0.0
!
!    m = 3.0
!
!    print *, "# time,  pos,  vel"
!    print *, time, pos, vel
!
!    do while (pos .ge. 0.0)
!        call update(time, dt,pos,vel,m,pos,vel)
!        time = time + dt
!        print *, time, pos, vel
!    end do

end program polytrope

!subroutine update(time,dt,x0,v0,m,x,v)
!    use solver, only : RK4
!    implicit none
!    external :: my_derive
!    real, intent(in)  :: time, dt, x0, v0, m
!    real, intent(out) :: x,v
!    integer,parameter :: n = 2
!    real,dimension(n) :: yin, ynext

!    yin(1) = x0
!    yin(2) = v0
!    call rk4(2,m,yin,ynext,time,dt, my_derive)

!    x  = ynext(1)
!    v  = ynext(2)
!    return
!end subroutine update

!subroutine my_derive(n,m,x,yin,k)
!    implicit none
!    integer, intent(in)  :: n  ! number of ODEs
!    real, intent(in)     :: m
!    real, intent(in)     :: x   !
!    real,dimension(n),intent(in)  :: yin   ! y
!    real,dimension(n),intent(out) :: k     ! dydx
!    k(1)  = yin(2)  ! v
!    if (x == 0) then
!        k(2) = -yin(1)**m     ! a
!   else
!        k(2) = -yin(1)**m-2.0/x*yin(2) 
!    end if
!    return
!end subroutine my_derive

