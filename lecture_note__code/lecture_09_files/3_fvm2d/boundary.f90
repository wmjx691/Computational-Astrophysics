subroutine boundary_x(v)
    use Simulation_data
    implicit none
    real,dimension(istart-ibuf:iend+ibuf,jstart-jbuf:jend+jbuf), intent(inout)  ::  v
    integer  ::  i,j
    
    ! left boundary
    do j = jstart, jend
        do i = 1, ibuf
            v(istart-i,j) = v(iend-i+1,j)
        enddo
    enddo
    
    ! right boundary
    do j = jstart, jend
        do i = 1, ibuf
            v(iend+i,j) = v(istart+i-1,j)
        enddo
    enddo
    
end subroutine boundary_x


subroutine boundary_y(v)
    use Simulation_data
    implicit none
    real,dimension(istart-ibuf:iend+ibuf,jstart-jbuf:jend+jbuf), intent(inout)  ::  v
    integer  ::  i,j
    
    ! down boundary
    do j = 1, jbuf
        do i = istart, iend
            v(i,jstart-j) = v(i,jend-j+1)
        enddo
    enddo
    
    ! top boundary
    do j = 1, jbuf
        do i = istart, iend
            v(i,jend+j) = v(i,jstart+j-1)
        enddo
    enddo
    
end subroutine boundary_y