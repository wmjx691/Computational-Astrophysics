program output
    implicit none
    type Star
        character(len=10) :: name
        integer(kind=4) :: id
        double precision :: age
        double precision :: distance
        double precision :: magnitude
        double precision :: is_double
    end type Star

    character*40         :: file_name
    character* 4         :: id_str
    integer,parameter    :: nstars = 100

    type(Star), dimension(nstars) :: stars
    integer :: n

    file_name = "stars txt"
    ! create the file
    open(unit=1,file=trim(File_name))
    
    ! write the header
    write(1,11) "# ", "Name", "ID", "double", "Age", "Distance", "Magnitude"
    write(1,11) "# ", " ", " ", " ", "[yr]", "[pc]", "[]"

    do n = 1, nstars
    ! generate the data
        write(id_str, 10) n
        stars(n)%name       =   "Star "//id_str
        stars(n)%id         =   n
        stars(n)%age        =   sqrt(real(n)**3)
        stars(n)%distance   =   real(n)**2 - 1.
        stars(n)%magnitude  =   real(n)**2.5
        stars(n)%is_double  =   (mod(n,5) .eq. 3)
        
        ! mrite to the file
        write(1,12) stars(n)%name, stars(n)%id, stars(n)%is_double, stars(n)%age, stars(n)distance, stars(n)%magnitude

    end do
    close(1)

    10  format(14)
    11  format(a2, a10, a5, a7, 3a24)
    12  format(2x, a10, i5, l7, 3e24.14)
    

end program output