subroutine initial(U)
implicit none
include 'parameters'
integer          :: i
double precision :: pi
pi = 4.d0*datan(1.d0)

do i=1,n
        U(i)= dsin(2.d0*pi*dble(i)/dble(n))
enddo

end
