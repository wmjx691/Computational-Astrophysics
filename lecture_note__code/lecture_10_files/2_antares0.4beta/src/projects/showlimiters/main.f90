program showlimiters
implicit none
include 'parameters'
integer          :: i
double precision :: fslop
double precision :: slop

! setupt initial condition
call initial(U)
! map U to UL
do i =1,n
        UL(3*i-1) = U(i)
enddo
! calculate slop and UL from limiter
do i =2,n-1
        slop = fslop(U(i-1),U(i),U(i+1))
        UL(3*i-2) = U(i)-0.5*slop
        UL(3*i  ) = U(i)+0.5*slop
enddo
open(3,file='showla.tab')
open(4,file='showlb.tab')
do i = 2,n-1
        write(3,*)dble(i)-0.5,UL(3*i-2) ! left  state
        write(3,*)dble(i)+0.5,UL(3*i  ) ! right state
        write(4,*)dble(i),UL(3*i-1)
enddo
close(3)
close(4)

write(*,*)'Finished!  '

end program
