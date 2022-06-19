module output
    implicit none

    contains

    subroutine record(time)
        use constants
        use Simulation_data
        implicit none
        real               :: a, e, L=0.0, energy
        real, intent(in)   :: time
        integer            :: i, j, k
        character(7)       :: starstr
        character(58)      :: filename

        logical, save      :: first = .true.


        a = 0.5*(au+0.5349192658994514E+14)
        e = au/a-1



11  format(i3,a4)

        do i=1, N
            write(starstr,11) i, '.dat'
            do j=1,3
                if(starstr(j:j).eq.' ') starstr(j:j)='0'
            enddo

            filename = 'binary_'//starstr

            if (first) then
                open(120,file=filename,status='unknown')
                write(120, 29)  "# Tag", "Time [sec]", "Mass [g]", &
                                "posx [code]", "posy [code]", &
                                "velx [code]", "vely [code]", &
                                "accx [code]", "accy [code]", &
                                "L [outcode]", "E [outcode]"
            else
                open(120,file=filename,status='old',position='append')
            endif

            !角動量是r x p
            do k=1, N
                L=L+(stars(k)%mass)*(stars(k)%x)*(stars(k)%vy)-(stars(k)%mass)*(stars(k)%y)*(stars(k)%vx)
            enddo
            energy=(0.5*(stars(2)%mass)*(G*(stars(1)%mass)*(stars(2)%mass)/L)**2)*(e**2-1)

            write(120, 30) i, time, stars(i)%mass, &
                            stars(i)%x, stars(i)%y, &
                            stars(i)%vx, stars(i)%vy, &
                            stars(i)%ax, stars(i)%ay, &
                            L, energy
            close(120)

            L=0

        enddo

        first = .false.

29  format(a6, 10a24)
30  format(i6, 10e24.12)
    
        return
    end subroutine record

end module output
