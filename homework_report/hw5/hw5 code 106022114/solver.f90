module solver

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    implicit none
    contains 

        !--------------------------------------------------
        !
        ! Implment the bisection method to solve the func
        !
        !
        ! Inputs:   func
        !
        ! Outputs:
        !           xs   : solution of func
        !           error: relative error
        !
        !--------------------------------------------------


        subroutine bisection(func, xs, err)
        implicit none
        real, external    :: func    ! the function to solve
        real, intent(out) :: xs      ! solution
        real, intent(out) :: err     ! error
        real, save :: a = -4.0        ! bracking interval [a,b]
        real, save :: b = -3.0        ! bracking interval [a,b]
        real  :: fa, fx              ! f(a) and f(x)

        xs=(a+b)/2

        fa=func(a)
        fx=func(xs)
        if (fa*fx<0) then
            b=xs
        else if (fa*fx>0) then
            a=xs          
        end if
        err=abs((b-a)/xs)

        end subroutine bisection


        !--------------------------------------------------
        !
        ! Implment the false position method to solve my_func
        !
        !
        ! Inputs: func
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine false(func, xs, err)
        implicit none
        real, external    :: func ! the function to solve
        real, intent(out) :: xs   ! solution
        real, intent(out) :: err  ! error

        real, save :: a = -4.0     ! bracking interval [a,b]
        real, save :: b = -3.0     ! bracking interval [a,b]
        real :: fa, fb            ! f(a) and f(b)
        real :: fx                ! f(xs)

        fa=func(a)
        fb=func(b)

        xs=a-fa*(b-a)/(fb-fa)
        fx=func(xs)

        if (fa*fx<0) then
            b=xs
        else if (fa*fx>0) then
            a=xs          
        end if
        err=abs(fx)
        

        end subroutine false


        !--------------------------------------------------
        !
        ! Implment the Newton's method to solve my_func
        !
        !
        ! Inputs: func, dfunc
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        ! subroutine newton(func, dfunc, xs, err)
        ! implicit none
        ! real, external    :: func  ! the function to solve
        ! real, external    :: dfunc ! the first derivative of the function to solve
        ! real, intent(out) :: xs    ! solution
        ! real, intent(out) :: err   ! error
        ! real, save :: x = 0.0      ! trial value
        ! real :: fx
        ! real :: dfdx

        ! fx=func(x)
        ! dfdx=dfunc(x)

        ! xs=x-fx/dfdx
        ! x=xs
        
        ! err=abs(fx)
        ! end subroutine newton

        subroutine solve_lower_triangular_matrix(N,L,b,x)
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: L  ! lower triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs              

            integer :: i,j

            bs=b
            do j=1,N
                if ( L(j,j) == 0 ) then
                    stop
                endif
                
                x(j) = bs(j)/L(j,j)
                
                do i=j+1,N
                    bs(i)=bs(i)-L(i,j)*x(j)
                enddo               
            enddo

            return
        end subroutine solve_lower_triangular_matrix

        subroutine solve_upper_triangular_matrix(N,U,b,x)
            implicit none

            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: U  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution
            real, dimension(N)  :: bs

            integer :: i,j

            bs=b
            do j=N,1,-1
                if ( U(j,j) == 0 ) then
                    stop
                endif
                
                x(j) = bs(j)/U(j,j)
                
                do i=1,j-1
                    bs(i)=bs(i)-U(i,j)*x(j)
                enddo               
            enddo

            return
        end subroutine solve_upper_triangular_matrix

        subroutine LU_decomposition(N,A,L,U)
            implicit none
            
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A    ! matrix
            real, dimension(N,N), intent(out)    :: L,U  ! matrix
            real, dimension(N,N) :: M, As

            integer :: i,j,k

            forall (j=1:N, i=1:N)
                L(i,j)=merge(1.0,0.0,i==j)
                U(i,j)=0.0
            end forall

            As=A
            do k=1,N-1
                if ( As(k,k) == 0 ) then
                    stop
                endif

                do i=k+1,N
                    M(i,k)=As(i,k)/As(k,k)
                enddo 
                
                do j=k+1,N
                    do i=k+1,N
                        As(i,j)=As(i,j)- M(i,k)*As(k,j)
                    enddo 
                enddo               
            enddo

            do i=1,N
                L(i, :i-1) = M(i, :i-1)
                U(i,i:   ) = As(i,i:  )
            enddo

        end subroutine LU_decomposition

        subroutine solve_lu(N,A,b,x)
            ! solve: A x = b
            integer, intent(in)  :: N
            real, dimension(N,N), intent(in)  :: A  ! upper triangle
            real, dimension(N)  , intent(in)  :: b  ! vector
            real, dimension(N)  , intent(out) :: x  ! solution

            real, dimension(N,N) :: L, U
            real, dimension(N)   :: y

            call LU_decomposition(N,A,L,U)
            call solve_lower_triangular_matrix(N,L,b,y)
            call solve_upper_triangular_matrix(N,U,y,x)


        end subroutine solve_lu

        subroutine newton_n(count, xs)
            implicit none
            integer, intent(in) :: count               ! count
            real, dimension(2), intent(out)  :: xs  ! solution

            real, dimension(2)   :: F   ! F
            real, dimension(2,2) :: J   ! J
            real, save, dimension(2)  :: x=(/1,0/)  ! trial value
            real, dimension(2)        :: s=(/1,0/)
            integer :: k=0
            

            do while (k<count)
                F(1) = x(1)+2*x(2)-2
                F(2) = x(1)**2+4*x(2)**2-4
            
                J(1,1) = 1
                J(1,2) = 2
                J(2,1) = 2*x(1)
                J(2,2) = 8*x(2)

                call solve_lu(2,J,-F,s)

                print *,xs
                xs=x+s
                x=xs
                k=k+1
            enddo
            
        end subroutine newton_n


        !--------------------------------------------------
        !
        ! Implment the Secant method to solve my_func
        !
        !
        ! Inputs: None
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        
        subroutine secant(func, xs, err)
            implicit none
            real, external    :: func  ! the function to solve
            real, intent(out) :: xs    ! solution
            real, intent(out) :: err   ! error

            real,save :: x0 = -3.0  ! initial guess
            real,save :: x1 = -4.0  ! initial guess

            real :: fx0, fx1, fxk  ! f(x0) and f(x1)

            fx0=func(x0)
            fx1=func(x1)

            fxk=((fx1-fx0)/(x1-x0))
            x0=x1

            x1=x1-fx1/fxk
            xs=x1

            err=abs(fx0)
        end subroutine secant

        subroutine mc_integral(func, a, b, N, answer)
            implicit none
            real,    external    :: func  ! the function to solve
            real,    intent(in)  :: a    ! lower
            real,    intent(in)  :: b    ! upper
            integer, intent(in)  :: N    ! number of calculation
            real,    intent(out) :: answer    ! solution

            real,    dimension(N):: xx,yy
            real,    dimension(N):: samples_x, samples_f
            integer              :: i,k
            real                 :: count,rangex=2.0, rangey=1.0

            do i=1,N
                xx(i)=(rangex/N)*i
                yy(i) = func(xx(i))
            enddo

            call random_number(samples_x) 
            samples_x = rangex*samples_x
            call random_number(samples_f)
            samples_f =  rangey*samples_f
            
            count=0.0
            
            do k=1,N
                if (yy(k)>=samples_f(k)) then 
                    count=count+1
                endif
            enddo
                    
            answer=(count/N)*(rangex*rangey)       
        
        end subroutine mc_integral

        subroutine integral(my_func, a, b, answer)
            implicit none
            real,    external    :: my_func  ! the function to solve
            real,    intent(in)  :: a    ! lower
            real,    intent(in)  :: b    ! upper
            real,    intent(out) :: answer    ! solution
        
            answer=(b-a)/6*(my_func(a)+4*my_func((a+b)/2)+my_func(b))

        end subroutine integral


end module solver
