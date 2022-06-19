!
! A simple program to calculate pi
!
program pi

    implicit none

    real    :: x, dx, h, dA, area
    integer :: i, N
    real    :: my_func

    print *, "Type the number of partition (N) you want to make"
    read *, N

    !-- initialize values
    dx   = (1.0 - (-1.0))/ N
    area = 0.0

    !-- calculate each rectangle and sum the area
    do i = 1, N
       x  = -1.0 + dx*i
       h  = my_func(x)
       dA = dx*h
       area = area + dA
    enddo

    !-- print out the result
    print *, "PI = ", 2.*area


end program pi

real function my_func(x)
    ! the function return the y values of a half circle with radius = 1
    real :: x
    my_func = sqrt(1.0-x**2)
    return
end function
