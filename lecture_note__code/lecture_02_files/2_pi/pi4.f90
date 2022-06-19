!
! A simple program to calculate pi
!
program calculate_pi

    implicit none

    real    :: area, error
    integer :: i, n
    real,parameter :: pi = 4.0*atan(1.0)
    real    :: t1, t2

    integer, parameter :: MAXN = 6
    integer, dimension(MAXN) ::  numbers

    !--- setup numbers
    numbers = (/10, 100, 1000, 10000, 100000, 1000000/)

    !--- open a file to store the error
    open(unit=11, file="error.dat")
    write(11,100) "# ", "N", "Error", "Compute time [s]"

    !--- measure the compute time
    call cpu_time(t1)
    ! do integral
    do i=1, MAXN
        n = numbers(i)
        call calculate_function_integral(n, area)
        print *, "N = ",n, " PI = ", 2.*area
        error = abs(2.*area - pi)/pi 

        call cpu_time(t2)
        write(11,200) numbers(i), error, (t2-t1)
        t1 = t2
    enddo

100 format(a2, a8, 2a24)
200 format(2x, i8, 2e24.14)

    close(11)

end program calculate_pi

subroutine calculate_function_integral(N, A)
    implicit none
    integer, intent(in)  :: N 
    real,    intent(out) :: A ! area of the function

    integer :: i
    real    :: my_func
    real    :: x, h, dx, dA

    !-- initialize values
    dx   = (1.0 - (-1.0))/ N
    A    = 0

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
