        program binary_out
        implicit none
        integer n,i,j,k,nend
        parameter (nend=32)
        real*8 U,R,head
        dimension  U(nend,nend,nend,8)
        dimension  R(196,196,1,8)
        
        do n =1,8
                do k=1,nend
                do j=1,nend
                do i=1,nend
                        u(i,j,k,n) = dble(i+j+k)*dble(n)
                enddo
                enddo
                enddo
        enddo
        open(60,file='test.bin',status='unknown',form='unformatted')
        write(60) nend,nend,nend,8
                do n=1,8
                do k=1,nend
                do j=1,nend
                        write(60) (U(i,j,k,n),i=1,nend)
                enddo
                enddo
                enddo

        close(60)

        open(60,file='dat0096.bin',status='unknown',form='unformatted')
        read(60) head
        write(*,*)head

                do n=1,8
                do k=1,1
                do j=1,196
                        read(60) (R(i,j,k,n),i=1,196)
                enddo
                enddo
                enddo
        close(60)
        !write(*,*) U-R
        !write(*,*) U
        write(*,*) R
        end program
