program data_type
    implicit none
    logical         ::  is_person
    character*40    ::  name
    complex         ::  complex_var
    real, parameter ::  msun = 1.989e33

    num_images  =   "kevin"
    is_person   =   .true
    complex_var =   (3.5)

    print *, name
    print *, is_person
    print *, complex_var
    print *, msun

end program data_type