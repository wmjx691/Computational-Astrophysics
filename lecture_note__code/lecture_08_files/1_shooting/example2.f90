program example2
    implicit none

    real :: dt, time, vel, pos

    real :: error

    real, save :: try, a, b
    integer :: iter, iter_max=100

    real    :: tol     = 1.e-8
    real    :: bound_pos   
    logical :: correct = .false.

    !上下界
    a    = -1.5
    b    = 2.0

    !期望的邊界條件
    bound_pos = 1.0
    iter      = 0

    do while (.not. correct)

        time = 0.0
        dt   = 0.01

        try  = 0.5*(a+b)
        pos  = 1.0
        vel  = try

        print *, "Iter = ", iter,"Trying y2 = ", try
        
        
        do while (time .le. 1.0)
            call update(time, dt, pos, vel, pos, vel)
            time = time + dt
        end do

        error = abs((pos-bound_pos)/bound_pos)
        if (error .le. tol) then
            correct = .true.
        endif

        if (pos .gt. bound_pos) then
            b = try
        else
            a = try
        endif

        iter = iter+1
        if (iter .ge. iter_max) then
            print *, "Error: Bisection method reached the max iteration"
            stop
        endif
    enddo
    
end program example2


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
    k(1)  = yin(2)  ! v
    k(2)  = 6.0*x     ! a
    return
end subroutine my_derive

