module IO
    implicit none

    contains

    subroutine record(time)
        use constants
        use Simulation_data
        implicit none
        real, intent(in)   :: time
        integer            :: i,j
        character(7)       :: starstr
        character(58)      :: filename

        logical, save      :: first = .true.

11  format(i3,a4)

        print *, "Output: time = ",time/yr, " [yr]"
        do i=1, N
            write(starstr,11) i, '.dat'
            do j=1,3
                if(starstr(j:j).eq.' ') starstr(j:j)='0'
            enddo

            filename = 'nbody_'//starstr

            if (first) then
                open(100,file=filename,status='unknown')
                write(100, 29)  "# Tag", "Time [sec]", "Mass [g]", &
                                "posx [code]", "posy [code]", &
                                "velx [code]", "vely [code]", &
                                "accx [code]", "accy [code]"
            else
                open(100,file=filename,status='old',position='append')
            endif

            write(100, 30) i, time, stars(i)%mass, &
                            stars(i)%x, stars(i)%y, &
                            stars(i)%vx, stars(i)%vy, &
                            stars(i)%ax, stars(i)%ay

            close(100)
        enddo

        first = .false.

29  format(a6, 8a24)
30  format(i6, 8e24.12)
    
        return
    end subroutine record


    subroutine read_model(names,masses,radii,angles)
        use Simulation_data
        implicit none
        integer :: i
        character*8,dimension(N),intent(out) :: names
        real, dimension(N), intent(out) :: masses, radii, angles

        print *, "reading model ...."
        open(unit=2,file='model.txt',status='old')
20  format(a8,3e10.3)

        do i = 1, N
            read(2,20) names(i), masses(i), radii(i), angles(i)
        enddo
        return
    end subroutine read_model
end module IO
