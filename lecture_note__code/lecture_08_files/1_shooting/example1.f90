program example1
    implicit none

    real :: dt, time, vel, pos

    dt   = 0.01  
    time = 0.0  

    pos  = 1.0
    vel  = -1.0

    print *, "# time,  pos,  vel"
    print *, time, pos, vel

    do while (time .le. 1.0)
        call update(time, dt,pos,vel,pos,vel)
        time = time + dt
        print *, time, pos, vel
    end do
end program example1


subroutine update(time,dt,x0,v0,x,v)
    use solver
    implicit none
    external :: my_derive
    real, intent(in)  :: time, dt, x0, v0
    real, intent(out) :: x,v
    integer,parameter :: n = 2
    real,dimension(n) :: yin, ynext

    yin(1) = x0
    yin(2) = v0
    call rk4(2, yin, ynext, time, dt, my_derive)

    x  = ynext(1)
    v  = ynext(2)
    return
end subroutine update

subroutine my_derive(n,x,yin,k)
    implicit none
    integer, intent(in)  :: n   ! number of ODEs
    real, intent(in)     :: x   !
    real,dimension(n),intent(in)  :: yin   ! y
    real,dimension(n),intent(out) :: k     ! dydx
    k(1)  = yin(2)  ! v
    k(2)  = 6.0*x     ! a
    return
end subroutine my_derive