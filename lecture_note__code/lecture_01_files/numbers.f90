program data_type
    implicit none
    integer*2 :: short
    integer*4 :: long
    integer*8 :: verylong
    integer*16 :: veryverylong

    real*4 :: real =1.0
    real*8 :: double =1.0
    real*16 :: quad =1.0

    print *, "Integers"

    print *,huge(short)
    print *,huge(long)
    print *,huge(verylong)
    print *,huge(veryverylong)

    print *, "Real numbers"
    print *,real
    print *,double
    print *,quad

end program data_type
