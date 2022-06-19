module IO
    implicit none
    contains
        subroutine output(n,time)
            use Simulation_data
            implicit none
            integer, intent(in) :: n
            real, intent(in)    :: time

            integer      :: i,j
            character(7) :: cycstr
            character(58)      :: filename

            real         :: ua ! analytical solution
            real         :: pa, pi, pe
            
            write(cycstr,10) n,'.d'
10  format(i5,a2)

            ! fill remaing zeros
            do i=1,5
                if(cycstr(i:i).eq.' ') cycstr(i:i)='0'
            enddo

            filename = 'advection_'//cycstr
            open(100,file=filename,status='unknown')

            ! write header
            ! write(100,29) "# i", "x", "U(x)", "UA(x)"

            do j=jend,jstart,-1
                do i=istart,iend
                    write(100,30) x(i), y(j), u(i,j)
                enddo
                write(100,*)" "
            enddo

29 format(a6,3a24)
30 format(3e24.12)

            close(100)
            return
        end subroutine output

end module IO
