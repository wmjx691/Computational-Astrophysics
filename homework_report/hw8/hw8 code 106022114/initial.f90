subroutine initial()
    use Simulation_data
    implicit none
    integer  :: i

    dx = 0.05   ! each point's x location

    !inital x array
    do i = istart, iend
        x(i)=xmin + i*dx      ! x will at cell center
    enddo

    !inital u
    do i = istart, iend
        if ((x(i) .ge. 0.0) .and. (x(i) .le. 0.5)) then
            u(i) = 2.0*x(i)
        else
            u(i) = 2.0-2.0*x(i)
        endif
    enddo

    return
end subroutine initial
