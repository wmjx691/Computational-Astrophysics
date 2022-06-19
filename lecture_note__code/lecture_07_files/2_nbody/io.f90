module IO
    implicit none

    contains

    subroutine record(time)
        use constants
        use Simulation_data
        implicit none
        real, intent(in)   :: time
        integer            :: i=1,j
        character(7)       :: starstr
        character(58)      :: filename

        logical, save      :: first = .true.

11  format(i3,a4)

        print *, "Output: time = ",time/yr, " [yr]"

            write(starstr,11) i, '.dat'
            do j=1,3
                if(starstr(j:j).eq.' ') starstr(j:j)='0'
            enddo

            filename = 'nbody_'//starstr

            if (first) then
                open(100,file=filename,status='unknown')
                write(100, 29)  "# Tag", "Time [sec]", &
                                "posx [code]", "posy [code]", &
                                "r [code]", &
                                "Angular momentum [code]", "energy [code]",&
                                "Jin [code]"
            else
                open(100,file=filename,status='old',position='append')
            endif

            write(100, 30) i, time, &
                            stars%r*cos(stars%theta), &
                            stars%r*sin(stars%theta), &
                            stars%r, &
                            stars%l , stars%energy, &
                            stars%Jin

            close(100)

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
