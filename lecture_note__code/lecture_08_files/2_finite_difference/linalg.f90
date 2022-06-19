module linalg
    ! copy your linalg.f90 from lecture 05
    implicit none
    contains 

    subroutine mat_print(amsg,a)
        character(*), intent(in) :: amsg
        class    (*), intent(in) :: a(:,:)
        integer                  :: i
        print*,' '
        print*,amsg
        do i=1,size(a,1)
            select type (a)
                type is (real(8)) ; print'(100f8.3)',a(i,:)
                type is (integer) ; print'(100i8  )',a(i,:)
            end select
        end do
        print*,' '
    end subroutine

    subroutine solve_lower_triangular_matrix(N,L,b,x)
        implicit none

        integer, intent(in)  :: N
        real*16, dimension(N,N), intent(in)  :: L  ! lower triangle
        real*16, dimension(N)  , intent(in)  :: b  ! vector
        real*16, dimension(N)  , intent(out) :: x  ! solution
        real*16, dimension(N)  :: bs              

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
        real*16, dimension(N,N), intent(in)  :: U  ! upper triangle
        real*16, dimension(N)  , intent(in)  :: b  ! vector
        real*16, dimension(N)  , intent(out) :: x  ! solution
        real*16, dimension(N)  :: bs

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
        real*16, dimension(N,N), intent(in)  :: A    ! matrix
        real*16, dimension(N,N), intent(out)    :: L,U  ! matrix
        real*16, dimension(N,N) :: M, As

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

    subroutine LU_decomposition_partial_pivoting(N,A,L,U,P)
        implicit none
        
        integer, intent(in)  :: N                      ! Ths sisze of the matrix
        real, dimension(N,N), intent(in)  :: A         ! matrix (N x N)
        real, dimension(N,N), intent(out)    :: L,U,P  ! matrix (N x N)
        real, dimension(N,N) :: M, As

        integer :: i,j,k,ip

    end subroutine LU_decomposition_partial_pivoting 

    subroutine swap_row(N,A,r1,r2)
        implicit none
        integer, intent(in) :: N
        integer, intent(in) :: r1, r2  ! swap row r1 and r2
        real, dimension(N,N), intent(inout)  :: A  
        real, dimension(N) :: row
    end subroutine swap_row

    subroutine solve_lu(N,A,b,x)
        ! solve: A x = b
        integer, intent(in)  :: N
        real*16, dimension(N,N), intent(in)  :: A  ! upper triangle
        real*16, dimension(N)  , intent(in)  :: b  ! vector
        real*16, dimension(N)  , intent(out) :: x  ! solution

        real*16, dimension(N,N) :: L, U
        real*16, dimension(N)   :: y

        call LU_decomposition(N,A,L,U)
        call solve_lower_triangular_matrix(N,L,b,y)
        call solve_upper_triangular_matrix(N,U,y,x)


    end subroutine solve_lu

    subroutine cholesky_factorization(N,A,L)
        implicit none
        integer, intent(in)  :: N                   ! Ths sisze of the matrix
        real, dimension(N,N), intent(in)  :: A      ! matrix (N x N)
        real, dimension(N,N), intent(out) :: L      ! matrix (N x N)

        integer :: i,j,k
        real, dimension(N,N) :: As

        As=A
        do k=1,N
            As(k,k)=As(k,k)**0.5
            L(k,k)=As(k,k)
            do i=k+1,N
                As(i,k)=As(i,k)/As(k,k)
                L(i,k)=As(i,k)
            enddo

            do j=k+1,N
                do i=k+1,N
                    As(i,j)=As(i,j)-As(i,k)*As(j,k)
                enddo
            enddo
        enddo        
    end subroutine cholesky_factorization

end module linalg
