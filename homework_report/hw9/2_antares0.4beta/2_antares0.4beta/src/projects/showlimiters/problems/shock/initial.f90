subroutine initial(U)
implicit none
include 'parameters'
integer          :: i
double precision :: pi
pi = 4.d0*datan(1.d0)

do i=1,n/2
        U(i)= dsin(1.d0*pi*dble(i)/dble(n))
enddo
do i=n/2,n
        U(i)= dcos(1.d0*pi*dble(i)/dble(n))
enddo

end
