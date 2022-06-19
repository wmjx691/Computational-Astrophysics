module Solver 
    ! copy your ODE solver from lecture 07
    implicit none
    contains
    subroutine rk4(n, yin, ynext, x, h, derive)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: x, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: derive

        integer :: i
        real              :: h2
        real,dimension(n) :: k1,k2,k3,k4
        real,dimension(n) :: y2,y3,y4

        call derive(n, x, yin, k1)
        do i=1,n
            y2(i)=yin(i)+h*k1(i)/3
        enddo
        call derive(n, x, y2, k2)
        do i=1,n
            y3(i)=y2(i)+h*k2(i)*2/3
        enddo
        call derive(n, x, y3, k3)
        do i=1,n
            y4(i)=y3(i)+h*k3(i)*2/3
        enddo
        call derive(n, x, y4, k4)
        do i=1,n
            ynext(i)=y4(i)+h*k4(i)/3
        enddo
        do i=1,n
            ynext(i)=0.5*(yin(i)+ynext(i))
        enddo
    end subroutine rk4
end module Solver
