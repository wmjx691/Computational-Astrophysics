subroutine initial()
    use Simulation_data
    implicit none
    integer  :: i

    dx = (xmax - xmin)/real(imax)

    !inital x array
    do i = istart-ibuf, iend+ibuf
        x(i)=xmin + (i-0.5)*dx
    enddo

    !inital u
    do i = istart, iend
        if ((x(i) .ge. 0.1) .and. (x(i) .le. 0.2)) then
            u(i) = 1.0
        else
            u(i) = 0.01
        endif
        !u(i)= exp(-((x(i)-0.1)*10.0)**2)
    enddo

    return
end subroutine initial
