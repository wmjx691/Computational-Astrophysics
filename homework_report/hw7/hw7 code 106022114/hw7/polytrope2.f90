!---------------------------------------------------
! The main program
!
program polytrope2
    implicit none
    integer, parameter  :: NMAX = 100    ! Max iteration number
    real, parameter     :: eps  = 1.e-6 ! tolerance

    integer  :: N              ! number of iteration
    real     :: error          ! error
    real     :: xs             ! solution
    real     :: dt, time, pos, m, vel, M_

    real, external :: func  ! the function to solve

    N     = 1
    error = 1e99
    do while (error > eps)
        call bisection(func, xs, error)
        !print *, "N = ",N, " solution is ", xs, " error = ", error
        N = N+1
        if (N .gt. NMAX) then
            print *, "The problem is not converged within ", N, " iterations."
            stop
        endif
    enddo
    
    dt   = 1000.0
    time = 1.0
    pos  = xs
    m = 1.0
    vel  = 0.0
    do while (time .le. 1200000.0)
        time = time + dt
        M_   = M_ + 4.0 * 3.1415926535 * time**2 * dt * pos
        call update(time, dt, pos, vel, m, pos, vel, M_)
        call record(time, pos, vel)
    end do

end program polytrope2

subroutine update(time,dt,x0,v0,m,x,v, M_)
    use solver, only : RK4
    implicit none
    external :: my_derive
    real, intent(in)  :: time, dt, x0, v0, m, M_
    real, intent(out) :: x,v
    integer,parameter :: n = 2
    real,dimension(n) :: yin, ynext

    yin(1) = x0
    yin(2) = v0
    call rk4(2,m,yin,ynext,time,dt, my_derive, M_)

    x  = ynext(1)
    v  = ynext(2)


    return
end subroutine update

subroutine my_derive(n,m,x,yin,k,M_)
    implicit none
    integer, intent(in)  :: n  ! number of ODEs
    real, intent(in)     :: m
    real, intent(in)     :: x   !
    real, intent(in)     :: M_
    real,dimension(n),intent(in)  :: yin   ! y
    real,dimension(n),intent(out) :: k     ! dydx
    
    real,parameter    ::  G = 6.67428e-8
    real,parameter    ::  pi= 3.1415926535
    real,parameter    ::  K_= 1.455e5
    real,parameter    ::  c = 2.99792458e10
    
    real  ::  ma
    
    ma = 4.0*pi*x**2*yin(1)

    k(1) =  yin(2)
            !-(G*M_/(2.0*K_)) &
            !*(M_*x**-2+(K_/c**2)*M_*yin(1)*x**-2+(4.0*pi*K_/c**2)*x*yin(1)**2+(4.0*pi*K_**2/c**4)*x*yin(1)**3)&
            !*(1.0-(2.0*G*M_/x*c**2))**(-1)
    
            !-(G*M_/(2.0*K_*x**2))*(1.0+K_*yin(1)/c**2) &
            !*(1.0+4.0*pi*x**3*K_*yin(1)**2/M_*c**2)*(1.0-(2.0*G*M_/x*c**2))**(-1)
   
    k(2) =  -(G*M_/(2.0*K_)) &
            *((ma*(x**-2)-2.0*M_*(x**-3))+(K_/c**2)*M_*yin(1)*x**-2+(4.0*pi*K_/c**2)*x*yin(1)**2+(4.0*pi*K_**2/c**4)*x*yin(1)**3 &
                        
            +(M_*x**-2+(K_/c**2)*(ma*yin(1)*x**-2+M_*yin(2)*x**-2-2.0*M_*yin(1)*x**-3)&
            +(4.0*pi*K_/c**2)*x*yin(1)**2+(4.0*pi*K_**2/c**4)*x*yin(1)**3)&
                        
            +(M_*x**-2+(K_/c**2)*M_*yin(1)*x**-2+(4.0*pi*K_/c**2)*(yin(1)**2+ 2.0*x*yin(1)*yin(2))&
            +(4.0*pi*K_**2/c**4)*x*yin(1)**3)&
                                    
            +(M_*x**-2+(K_/c**2)*M_*yin(1)*x**-2+(4.0*pi*K_/c**2)*x*yin(1)**2+(4.0*pi*K_**2/c**4)&
            *(yin(1)**3+3.0*x*yin(1)**2*yin(2))))&
            
            *(1.0-(2.0*G*M_/x*c**2))**(-1) &
            
            -(G*M_/(2.0*K_)) &
            *(M_*x**-2+(K_/c**2)*M_*yin(1)*x**-2+(4.0*pi*K_/c**2)*x*yin(1)**2+(4.0*pi*K_**2/c**4)*x*yin(1)**3)&
            *((1.0-(2.0*G*M_/x*c**2))**(-2))*(2.0*G/c**2)*(ma*(x**-1)-M_*(x**(-2)))
           
    return
end subroutine my_derive

subroutine bisection(func, xs, err)
    implicit none
    real, external    :: func    ! the function to solve
    real, intent(out) :: xs      ! solution
    real, intent(out) :: err     ! error
    real, save :: a = 0.0       ! bracking interval [a,b]
    real, save :: b = 1.0e16   ! bracking interval [a,b]
    real  :: fa, fx              ! f(a) and f(x)

    xs = (a+b)*0.5
    fa = func(a,0.0)
    fx = func(xs,1.0)
    if (sign(1.0,fa) .eq. sign(1.0,fx)) then
        b = xs
    else
        a = xs
    endif
    err=abs(fx)
end subroutine bisection

real function func(x,i)
    use solver
    implicit none
    real, intent(in)     :: x,i
    real :: dt, time, vel, pos, m, M_
    dt   = 1.0
    time = 1.0
    pos  = x
    m = 1.0
    vel  = 0.0
    M_   = 0.0
    do while (time .le. 1200000.0)
        time = time + dt
        M_   = M_ + 4.0 * 3.1415926535 * time**2 * dt * pos
        if (i==1.0) then
           !print *,time,pos,vel
        endif
        call update(time, dt,pos,vel,m,pos,vel,M_)
    end do
    
    func = pos-0.0
    
    if (i == 1.0) then
        print *,"p=",pos,"p'=",vel,"error=",func
    endif
        
    
    
    return 
end function func

subroutine record(time,p,p_)
      implicit none
      real, intent(in)   :: time,p,p_
      integer            :: i=1,j
      character(7)       :: starstr
      character(58)      :: filename

      logical, save      :: first = .true.

11  format(i3,a4)

      !print *, "Output: radius = ",time, " [cm]"

          write(starstr,11) i, '.dat'
          do j=1,1
              if(starstr(j:j).eq.' ') starstr(j:j)='0'
          enddo

          filename = 'nbody_'//starstr

          if (first) then
              open(100,file=filename,status='unknown')
              write(100, 29)  "# Tag", "Radius [cm]", &
                              "p [code]", "p_ [code]"
          else
              open(100,file=filename,status='old',position='append')
          endif

          write(100, 30) i, time, &
                          p, &
                          p_
          close(100)

      first = .false.

29  format(a6, 8a24)
30  format(i6, 8e24.12)
    
      return
end subroutine record
