! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Kuo-Chuan Pan 2020.04.14
!
! Problem:
!
!        Solving non-linear equations
!
!
program nonlinear

    use solver
    implicit none

    integer, parameter  :: N = 2
    real, dimension(N) :: x
    real, dimension(N) :: inx=(/10,0/)
    integer :: count

    count=10

    print *, "The progress of finding solution..."

    call newton_n(count, x)

    print *, "we start from inital x =",inx

    print *, "solution x = ",x

    ! integer, parameter  :: NMAX = 50    ! Max iteration number
    ! real, parameter     :: eps  = 1.e-6 ! tolerance

    ! integer  :: N              ! number of iteration
    ! real     :: error          ! error
    ! real     :: xs             ! solution
    ! character*40  :: fname     ! output filename
    ! real, external :: my_func  ! the function to solve
    ! real, external :: my_dfunc ! the d function to solve

!     ! output file name
!     !fname = "bisection.txt"
!     !fname = "false.txt"
!     fname = "newton1.txt"
!     !fname = "secant.txt"

!     open(unit=1,file=trim(fname))
!     ! write the header
!     write(1,11) "#", "N", "solution","Error"

!     ! The main iteration loop
!     N     = 1
!     error = 1e99

!     do while (abs(error) > eps)
!         !call bisection(my_func, xs, error)
!         !call false(my_func, xs, error)
!         call newton(my_func, my_dfunc, xs, error)
!         !call secant(my_func, xs, error)

!         print *, "N = ",N, " solution is ", xs, " error = ", error
!         write(1,12) N, xs, error

!         N = N+1
!         ! check if we have reached the maximum iteration number
!         if (N .gt. NMAX) then
!             print *, "The problem is not converged within ", N, " iterations."
!             stop
!         endif
!     enddo

!     close(1)

! 11  format(a2,a4,2a24)
! 12  format(2x,i4,2e24.14)

!     integer  :: N=1000     ! number of iteration
!     real     :: lower=0
!     real     :: upper
!     real     :: xs             ! solution
!     real     :: A, answer, error, realval=5.67037e-8

!     real, parameter     :: k  = 1.38e-23 ! tolerance
!     real, parameter     :: pi = 3.1415926535! tolerance
!     real, parameter     :: h  = 6.626e-34 ! tolerance
!     real, parameter     :: c  = 3.e8 ! tolerance

!     real, external :: my_func  ! the function to integral

!     upper=huge(upper)

!     call mc_integral(my_func, lower, upper, N, xs)

!     print *, "Number of N:"
!     print *, N
    
!     print *, "Result of the integral:"
!     print *, xs

!     A=(2.0*pi*k**4)/(h**3*c**2)

!     print *, "The coefficient of A:"
!     print *, A

!     answer=A*xs

!     print *, "The Stefan-Blotzmann constant is:"
!     print *, answer

!     error=(abs(answer-realval)/realval)*100

!     print *, "The error of 1e3 is %:"
!     print *, error

!     integer  :: N=1000     ! number of iteration
!     real     :: lower=0
!     real     :: upper
!     real     :: xs             ! solution
!     real     :: A, answer, error, realval=5.67037e-8

    ! real,    external    :: my_func
    ! real     :: dx=0.0000001    ! number of iteration
    ! real     :: lower=0.0
    ! real     :: upper,dA
    ! real     :: A=0.0         ! solution

    ! upper = lower+dx

    ! do while(upper<=2.0)
    !     call integral(my_func, lower, upper, dA)
    !     A      = dA+A
    !     lower  = dx+lower
    !     upper  = dx+upper
    ! enddo

    ! print *, "Result of the integral 0-2:"
    ! print *, A
    ! print *, "Result of the integral 2-infinite:"
    ! print *,0.03559
    ! A=A+0.03559
    ! print *, "Result of doing step integral 0-infinite:"
    ! print *, A


    ! print *, "Done!"
end program nonlinear 

!--------------------------------------------------
!
! my_function : The function to solve
!
!        f(x) = x^2 - 4 sin(x) = 0
!
!--------------------------------------------------
real function my_func(x)
        implicit none
        real, intent(in)  :: x
        real ::temp

        ! TODO:
        my_func=1/(x**4+x**2+1)
        return
end function my_func


! ! ---------------------------------------
! ! return f'(x) for Newton's method
! ! ---------------------------------------
! real function my_dfunc(x)
!     implicit none
!     real :: x

!     ! TODO:
!     my_dfunc=3*x**2+3*x-5.75

!     return
! end function my_dfunc



