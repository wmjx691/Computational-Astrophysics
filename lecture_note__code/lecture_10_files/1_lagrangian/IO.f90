module IO
    implicit none
    contains
        subroutine output(n,time)
            use Simulation_data
            implicit none
            integer, intent(in) :: n
            real, intent(in)    :: time
            real ,parameter     :: small = 1.e-99

            integer      :: i, ix
            character(7) :: cycstr
            character(58)      :: filename

            
            write(cycstr,10) n,'.d'
10  format(i5,a2)

            ! fill remaing zeros
            do i=1,5
                if(cycstr(i:i).eq.' ') cycstr(i:i)='0'
            enddo

            filename = 'hydro_'//cycstr
            open(100,file=filename,status='unknown')

            ! write header
            write(100,101) "# i", "fm", "r/r(nx)", "u", "rho", "p", "eps", "w"

            write(100,102)(ix,max(fm(ix),small),max(r(ix)/r(nx),small),max(u(ix),small), &
                max(rho(ix),small),max(p(ix),small),max(eps(ix),small),max(w(ix),small),ix=1,nx)

101 format(1x,a5,7a13)
102 format(1x,i5,1p7e12.4)

            close(100)
            return
        end subroutine output

end module IO
