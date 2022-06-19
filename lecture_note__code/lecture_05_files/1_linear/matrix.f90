!
!
! 
!
program linear
    use linalg
    implicit none
    integer, parameter  :: N = 2
    real*16,dimension(N,N) :: lower, upper, A, P, Ainv
    real*16,dimension(N) :: b
    real*16,dimension(N) :: x
    real,dimension(4,4) :: aa,ll,uu,pp
    real*16 :: r
    integer :: i,j

    r=1d-6

    A(1,1) =  1.0
    A(1,2) =  1.0+r
    ! A(1,3) =  0.0

    A(2,1) =  1.0-r
    A(2,2) =  1.0
    ! A(2,3) =  0.0

    ! A(3,1) =  1.0
    ! A(3,2) =  2.0
    ! A(3,3) =  2.0

    ! the vectore b
    b(1) =  1.0+r+r**2
    b(2) =  1.0
    ! b(3) =  3.0

    call LU_decomposition(N,A,lower,upper)
    call mat_print("A",A)
    call mat_print("L",lower)
    call mat_print("U",upper)

    print *, "vector   b = ",b

    call solve_lu(N,A,b,x)
    print *, "solution x = ",x

    print *, "error of x0(1) = ",(1-x(1))/1

    print *, "error of x1(epsilon) = ",(x(2)-r)/r

    ! call solve_upper_triangular_matrix(N,upper,b,x)
    ! print *, "solution x = ",x

    ! call mat_print("A:",A)
    ! call mat_print("b:",b)
    ! print *, "small parameter we set = ",e
    ! call cholesky_factorization(N,A,lower)
    ! call mat_print("Using cholesky_factorization:",lower)
    ! call solve_lower_triangular_matrix(N,lower,b,x)
    ! print *, "solution x = ",x


end program linear 

