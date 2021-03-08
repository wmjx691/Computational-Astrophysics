!
! A simple program to calculate pi
!
program pi

    implicit none

    real    :: x, dx, h, dA, area
    integer :: i, N

    print *, "Type the number of partition (N) you want to make"
    read *, N

    !-- initialize values
    !   dx = (1 - (-1))/ N
    !
    dx   = (1.0 - (-1.0))/ N
    area = 0.0

    !-- calculate each rectangle and sum the area
    do i = 1, N
       x  = -1.0 + dx*i
       h  = sqrt(1.0-x**2.0)
       dA = dx*h
       area = area + dA
    enddo

    !-- print out the result
    print *, "PI = ", 2.*area


end program pi


