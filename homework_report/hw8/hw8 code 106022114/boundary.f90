subroutine bounary(v)
    use Simulation_data
    implicit none
    real,dimension(istart-ibuf:iend+ibuf), intent(inout)  ::  v
    integer  ::  i
    
    ! left boudndary
    do i = 1, ibuf
        v(istart-ibuf) = 0.0
        v(istart) = 0.0
    enddo
    
    ! right boundary
    do i = 1, ibuf
        v(iend+ibuf) = 0.0
        v(iend) = 0.0
    enddo
    
end subroutine