program polytrope
    use Solver, only : rk4
    implicit none

    real :: dt, time, y1, y2

    dt   = 0.01  
    time = 0.01

    y1  = 1.0
    y2  = 0.00001

    print *, "# time,  y1,  y2"
    print *, time, y1, y2

    do while (time .le. 5.0)
        if (y1 < 0.0) then
            stop
        endif

        call update(time, dt, y1, y2, y1, y2)
        time = time + dt
        print *, time, y1, y2
    end do


end program polytrope

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
    call rk4(2,yin,ynext,time,dt, my_derive)

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
    k(1)  = yin(2)                         ! v
    k(2)  = -yin(1)**1.0-(2.0/x)*yin(2)      ! a
    return
end subroutine my_derive