!
!
! 
!
program linear
    use linalg
    implicit none
    integer, parameter  :: N = 4
    real*16,dimension(N,N) :: lower, upper, A, P, Ainv
    real*16,dimension(N) :: b
    real*16,dimension(N) :: x
    real*16 :: r,h
    integer :: i,j

    r = 1d-6
    h = 0.2

    A(1,1) =  -2.0
    A(1,2) =  1.0
    A(1,3) =  0.0
    A(1,4) =  0.0

    A(2,1) =  1.0
    A(2,2) =  -2.0
    A(2,3) =  1.0
    A(2,4) =  0.0

    A(3,2) =  1.0
    A(3,3) =  -2.0
    A(3,4) =  1.0
    A(3,1) =  0.0

    A(4,3) =  1.0
    A(4,4) =  -2.0
    A(4,1) =  0.0
    A(4,2) =  0.0

    ! the vectore b
    b(1) =  4.8*h**2-1.0
    b(2) =  3.6*h**2
    b(3) =  2.4*h**2
    b(4) =  1.2*h**2-1.0

    call LU_decomposition(N,A,lower,upper)
    call mat_print("A",A)
    call mat_print("L",lower)
    call mat_print("U",upper)

    print *, "vector   b = ",b

    call solve_lu(N,A,b,x)
    print *, "solution x = ",x


end program linear 

