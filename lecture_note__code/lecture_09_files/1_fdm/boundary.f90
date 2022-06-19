subroutine bounary(v)
    use Simulation_data
    implicit none
    real,dimension(istart-ibuf:iend+ibuf), intent(inout)  ::  v
    integer  ::  i
    
    ! left boudndary
    do i = 1, ibuf
        v(istart-i) = v(iend-i+1)
    enddo
    
    ! right boundary
    do i = 1, ibuf
        v(iend+i) = v(istart+i-1)
    enddo
    
end subroutine