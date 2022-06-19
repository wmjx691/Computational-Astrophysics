!
! A simple program to calculate pi
!
program pi

    implicit none

    real    :: area
    integer :: N

    print *, "Type the number of partition (N) you want to make"
    read *, N

    ! do integral
    call calculate_function_integral(N, area)

    !-- print out the result
    print *, "PI = ", 2.*area

end program pi

subroutine calculate_function_integral(N, A)
    implicit none
    integer, intent(in)  :: N
    real,    intent(out) :: A

    integer :: i
    real    :: my_func
    real    :: x, h, dx, dA

    !-- initialize values
    dx   = (1.0 - (-1.0))/ N
    A    = 0.0

    !-- calculate each rectangle and sum the area
    do i = 1, N
        x  = -1.0 + dx*i
        h  = my_func(x)
        dA = dx*h
        A = A + dA
    enddo

    return
end subroutine calculate_function_integral

real function my_func(x)
    ! the function return the y values of a half circle with radius = 1
    real :: x
    my_func = sqrt(1.0-x**2)
    return
end function
