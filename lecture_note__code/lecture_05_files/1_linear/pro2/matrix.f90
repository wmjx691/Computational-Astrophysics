!
!
! 
!
program linear
    use linalg
    implicit none
    integer, parameter  :: N = 3
    real,dimension(N,N) :: lower, upper, A, P, Ainv
    real,dimension(N) :: b
    real,dimension(N) :: x
    real,dimension(4,4) :: aa,ll,uu,pp
    integer :: i,j

    A(1,1) =  4.0
    A(1,2) = -1.0
    A(1,3) = -1.0

    A(2,1) = -1.0
    A(2,2) =  4.0
    A(2,3) = -1.0

    A(3,1) = -1.0
    A(3,2) = -1.0
    A(3,3) =  4.0

    ! the vectore b
    b(1) =  2.0
    b(2) =  8.0
    b(3) =  10.0

    ! call LU_decomposition(N,A,lower,upper)
    ! call mat_print("A",A)
    ! call mat_print("L",lower)
    ! call mat_print("U",upper)

    ! print *, "vector   b = ",b

    ! call solve_lower_triangular_matrix(N,lower,b,x)
    ! print *, "solution x = ",x

    ! call solve_upper_triangular_matrix(N,upper,b,x)
    ! print *, "solution x = ",x

    call cholesky_factorization(N,A,lower)
    call mat_print("cholesky_factorization of pro4 in writing assignment",lower)
    
end program linear 

