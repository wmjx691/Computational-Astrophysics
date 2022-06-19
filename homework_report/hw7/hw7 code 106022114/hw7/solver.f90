module solver 

    implicit none
    contains

    subroutine euler(n, yin, ynext, x, h, derive)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: x, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: derive
        integer            :: i
        real, dimension(n) :: k1
        call derive(n,x,yin,k1)
        do i = 1,n
            ynext(i)=yin(i)+h*k1(i)    
        end do        
    end subroutine euler

    subroutine rk2(n, yin, ynext, x, h, derive)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: x, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: derive
        integer            :: i
        real, dimension(n) :: k1, k2
        real,dimension(n)  :: y2
        call derive(n,x,yin,k1)
        do i = 1,n
            ynext(i)=yin(i)+h*k1(i)
        end do

        call derive(n,x+h,ynext,k2)
        do i = 1,n
            ynext(i)=yin(i)+h/2*(k1(i)+k2(i))
        end do

    end subroutine rk2

    subroutine rk4(n, m,  yin, ynext, x, h, derive, M_)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: x, h
        real, intent(in)    :: m
	real, intent(in)    :: M_
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: derive

        integer :: i
        real              :: h2
        real,dimension(n) :: k1,k2,k3,k4
        real,dimension(n) :: y2,y3,y4

        call derive(n,m,x,yin,k1,M_)
        call derive(n,m,x+h/2.0,yin+h/2.0*k1,k2,M_)
        call derive(n,m,x+h/2.0,yin+h/2.0*k2,k3,M_)
        call derive(n,m,x+h,yin+h*k3,k4,M_)
        !do i=i,n
        !    y2(i)=yin(i)+h/2.0*k1(i)
        !end do
        !call derive(n,x+h/2.0,y2,k2)
        !do i=i,n
        !    y3(i)=yin(i)+h/2.0*k2(i)
        !end do
        !call derive(n,x+h/2.0,y3,k3)
        !do i=i,n
        !    y4(i)=yin(i)+h*k3(i)
        !end do
        !call derive(n,x+h,y4,k4)
        ynext=yin+h/6.0*(k1+2.0*k2+2.0*k3+k4)
        !do i=i,n
        !    ynext(i)=yin(i)+h/6.0*(k1(i)+2.0*k2(i)+2.0*k3(i)+k4(i))
        !end do
    end subroutine rk4
end module solver
