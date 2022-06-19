!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        subroutine initial()
            use constants, only : au, msun, pi, G
            implicit none
            integer :: i
            real :: m1, m2, m3, force, force3
            real :: v3, period3, separation3

        end subroutine initial

        subroutine initial_from_model()
                use IO, only : read_model
                use Simulation_data

                implicit none
                integer :: i

                character*8, dimension(N) :: names
                real,dimension(N) :: masses, radii, angles

                call read_model(names,masses,radii,angles)
                call convert(masses, radii, angles)

        end subroutine initial_from_model

        !!-----------------------------------------------------------
        ! convert masses [g] and radii [AU] into Cartesian cooridnate
        !
        ! We assume a random initial theta of position and circular orbit
        ! based on Kepler's law
        !
        !!-----------------------------------------------------------
        subroutine convert(masses,radii, angles)

            use constants, only : au, msun, pi, G
            use Simulation_data
            implicit none
            real,dimension(N),intent(in) :: masses, radii,angles

            integer*4,parameter :: seed = 9527
            integer             :: i
            real                :: tht, r, v, acc, p
            real, parameter     :: small = 1.e-5
            !!

            call srand(seed)
            v = 10**5

            !! assume circular orbit
            do i=1, N

                r   = radii(i)*au              ! [cm]
                tht = rand()*pi*2.0            ! [rad]

                p = (((4.0*pi**2)/(G*(msun+stars(i)%mass)))*r**3)**0.5 !TODO

                ! TODO:
                stars(i)%mass = masses(i)
                stars(i)%x    = r*cos(tht)
                stars(i)%y    = r*sin(tht)

                if(r <= 0.0) then
                    v   = 0.0
                    acc = 0.0
                else
                    v   = 2.0*pi*r/p
                    acc = -G*msun/r**2
                endif

                stars(i)%vx   = -v*sin(tht)
                stars(i)%vy   = v*cos(tht)
                stars(i)%ax   = acc*cos(tht)
                stars(i)%ay   = acc*sin(tht)

            enddo

            !!
            stars2 = stars ! mainly for mass
            stars3 = stars
            stars4 = stars
        end subroutine convert


        subroutine update(time,dt)
            use constants
            implicit none
            real, intent(in)  :: time, dt
            integer :: i
            real    :: h

            !euler's method
            !h=dt
            !call update_euler(h,stars,stars)
            !call update_acc(dt,stars,stars)

            ! stars2 = stars ! mainly for mass
            ! stars3 = stars
            ! stars4 = stars

            !! RK2

            !TODO:


            !! RK4

            ! TODO:
            h = 0.5*dt
            call update_acc(dt,stars,stars)          !k1
            call update_euler(h,stars,stars,stars2)  !y2

            h = 0.5*dt
            call update_acc(dt,stars2,stars2)        !k2
            call update_euler(h,stars,stars2,stars3) !y3

            h = dt
            call update_acc(dt,stars3,stars3)        !k3
            call update_euler(h,stars,stars3,stars4) !y4
            call update_acc(h,stars4,stars4)         !k4
     


            do i=1, N
            stars(i)%x = stars(i)%x + &
                         dt*(stars(i)%vx + 2.0*stars2(i)%vx &
                                         + 2.0*stars3(i)%vx &
                                         + 1.0*stars4(i)%vx)/6.0

            stars(i)%y = stars(i)%y + &
                         dt*(stars(i)%vy + 2.0*stars2(i)%vy &
                                         + 2.0*stars3(i)%vy &
                                         + 1.0*stars4(i)%vy)/6.0
       
            stars(i)%vx = stars(i)%vx + &
                         dt*(stars(i)%ax + 2.0*stars2(i)%ax &
                                         + 2.0*stars3(i)%ax &
                                         + 1.0*stars4(i)%ax)/6.0
       
            stars(i)%vy = stars(i)%vy + &
                         dt*(stars(i)%ay + 2.0*stars2(i)%ay &
                                         + 2.0*stars3(i)%ay &
                                         + 1.0*stars4(i)%ay)/6.0
            enddo

            return
        end subroutine update

        subroutine update_euler(dt,s1,s2,s3)
            !!
            !! use s1 and s2 -> s3
            !!
            use constants
            implicit none
            real, intent(in) :: dt
            type(Star), dimension(N),intent(in)  :: s1
            type(Star), dimension(N),intent(in)  :: s2
            type(Star), dimension(N),intent(out) :: s3

            integer :: i,j
            real    :: x, y, rsq, fx, fy
            real    :: radius, force, angle

            ! update position and velocity
            do i=1, N

                ! TODO:
                s3(i)%x  = s1(i)%x  +s2(i)%vx *dt
                s3(i)%y  = s1(i)%y  +s2(i)%vy *dt
                s3(i)%vx = s1(i)%vx +s2(i)%ax *dt
                s3(i)%vy = s1(i)%vy +s2(i)%ay *dt
            enddo

        end subroutine update_euler


        subroutine update_acc(dt,s1,s2)
            !!
            !! use s1 -> s2
            !!
            use constants
            implicit none
            real, intent(in) :: dt
            type(Star), dimension(N),intent(in)  :: s1
            type(Star), dimension(N),intent(out) :: s2

            integer :: i,j
            real    :: x, y, rsq, fx, fy
            real    :: radius, force, angle
            ! update accelerations
            do i=1,N
                fx=0
                fy=0
                do j=1,N
                    if (i<j) then
                        x = s1(i)%x - s1(j)%x
                        y = s1(i)%y - s1(j)%y
                    else if (i>j) then
                        x = -(s1(i)%x - s1(j)%x)
                        y = -(s1(i)%y - s1(j)%y)
                    else if (i==j) then
                        cycle
                    endif

                    rsq = x*x + y*y
                    angle = atan(y,x)
                    force = G *s1(i)%mass *s1(j)%mass /rsq

                    if (i<j) then
                        fx=fx-force*cos(angle)
                        fy=fy-force*sin(angle)

                    else if (i>j) then
                        fx=fx+force*cos(angle)
                        fy=fy+force*sin(angle)
                    endif

                enddo
                s2(i)%ax = fx/s1(i)%mass
                s2(i)%ay = fy/s1(i)%mass
            enddo

        end subroutine update_acc

end module physics

