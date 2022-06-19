!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        ! subroutine initial()
        !     use constants, only : au, msun, pi, G
        !     implicit none
        !     integer :: i
        !     real :: m1, m2, m3, force, force3
        !     real :: v3, period3, separation3

        ! end subroutine initial

        subroutine initial_from_model(elipson,omegap,vc,vr,m)
                use IO, only : read_model
                use Simulation_data

                implicit none
                integer :: i
                real, intent(in) :: elipson
                real, intent(in) :: omegap
                real, intent(in) :: vc,vr
                integer, intent(in) :: m


                character*8, dimension(N) :: names
                real,dimension(N) :: masses, radii, angles

                ! call read_model(names,masses,radii,angles)
                call convert(elipson,omegap,vc,vr,m)

        end subroutine initial_from_model

        ! !!-----------------------------------------------------------
        ! ! convert masses [g] and radii [AU] into Cartesian cooridnate
        ! !
        ! ! We assume a random initial theta of position and circular orbit
        ! ! based on Kepler's law
        ! !
        ! !!-----------------------------------------------------------
        subroutine convert(elipson,omegap,vc,vr,m)

            use constants, only : au, msun, pi, G
            use Simulation_data
            implicit none
            ! real,dimension(N),intent(in) :: masses, radii,angles

            ! integer*4,parameter :: seed = 9527
            integer             :: i
            real, intent(in) :: elipson
            real, intent(in) :: omegap
            real, intent(in) :: vc,vr
            integer, intent(in) :: m
            ! real                :: tht, r, v, acc, p
            ! real, parameter     :: small = 1.e-5
            !!

            ! call srand(seed)
            ! v = 10**5

            !! assume circular orbit

            stars%r     = 1.0   ! [cm]
            stars%theta = 0.0
            ! tht =                ! [rad]

            ! p = (((4.0*pi**2)/(G*(msun+stars(i)%mass)))*r**3)**0.5 !TODO

            ! TODO:
            ! stars(i)%mass = masses(i)
            ! stars(i)%x    = r*cos(tht)
            ! stars(i)%y    = r*sin(tht)

            ! if(r <= 0.0) then
            !     v   = 0.0
            !     acc = 0.0
            ! else
            !     v   = 2.0*pi*r/p
            !     acc = -G*msun/r**2
            ! endif


            stars%r1     = vr
            stars%theta1 = vc/stars%r
            stars%r2     = 0.0
            stars%theta2 = 0.0

            stars%J      = 0.5-omegap*1.0
            stars%Jin    = 0.0



            !!
            stars2 = stars ! mainly for mass
            stars3 = stars
            stars4 = stars
        end subroutine convert


        subroutine update(time,elipson,omegap,vc,m,dt)
            use constants
            implicit none
            real, intent(in)  :: time, dt
            real, intent(in) :: elipson
            real, intent(in) :: omegap
            real, intent(in) :: vc
            integer, intent(in) :: m
            integer :: i
            real    :: h

            !euler's method
            ! h=dt
            ! call update_euler(h,stars,stars,stars)
            ! call update_acc(dt,elipson,omegap,vc,m,stars,stars)

            ! stars2 = stars ! mainly for mass
            ! stars3 = stars
            ! stars4 = stars

            !! RK2

            !TODO:


            !! RK4

            ! TODO:
            h = 0.5*dt
            call update_acc(dt,elipson,omegap,vc,m,stars,stars)          !k1
            call update_euler(h,stars,stars,stars2)                      !y2

            h = 0.5*dt
            call update_acc(dt,elipson,omegap,vc,m,stars2,stars2)        !k2
            call update_euler(h,stars,stars2,stars3)                     !y3

            h = dt
            call update_acc(dt,elipson,omegap,vc,m,stars3,stars3)        !k3
            call update_euler(h,stars,stars3,stars4)                     !y4
            call update_acc(h,elipson,omegap,vc,m,stars4,stars4)         !k4
     


            stars%r = stars%r + &
                            dt*(stars%r1 + 2.0*stars2%r1 &
                                         + 2.0*stars3%r1 &
                                         + 1.0*stars4%r1)/6.0

            stars%theta = stars%theta + &
                            dt*(stars%theta1 + 2.0*stars2%theta1 &
                                             + 2.0*stars3%theta1 &
                                             + 1.0*stars4%theta1)/6.0
       
            stars%r1 = stars%r1 + &
                            dt*(stars%r2 + 2.0*stars2%r2 &
                                         + 2.0*stars3%r2 &
                                         + 1.0*stars4%r2)/6.0
       
            stars%theta1 = stars%theta1 + &
                            dt*(stars%theta2 + 2.0*stars2%theta2 &
                                             + 2.0*stars3%theta2 &
                                             + 1.0*stars4%theta2)/6.0

            return
        end subroutine update

        subroutine update_euler(dt,s1,s2,s3)
            !!
            !! use s1 and s2 -> s3
            !!
            use constants
            implicit none
            real, intent(in) :: dt
            type(Star),intent(in)  :: s1
            type(Star),intent(in)  :: s2
            type(Star),intent(out) :: s3

            integer :: i
            ! real    :: x, y, rsq, fx, fy
            ! real    :: radius, force, angle

            ! update position and velocity
            ! TODO:
            s3%r      = s1%r      +s2%r1 *dt
            s3%theta  = s1%theta  +s2%theta1 *dt
            s3%r1     = s1%r1     +s2%r2 *dt
            s3%theta1 = s1%theta1 +s2%theta2 *dt

        end subroutine update_euler 

        subroutine update_acc(dt,elipson,omegap,vc,m,s1,s2)
            !!
            !! use s1 -> s2
            !!
            use constants
            implicit none
            real, intent(in) :: dt
            real, intent(in) :: elipson
            real, intent(in) :: omegap
            real, intent(in) :: vc
            integer, intent(in) :: m

            type(Star),intent(in)  :: s1
            type(Star),intent(out) :: s2

            ! integer :: i,j
            ! real    :: x, y, rsq, fx, fy
            real    :: angle,fr,ft,temp

            angle=m*(s1%theta-omegap*dt)*pi*2.0

            fr = -(vc**2/(s1%r))*(1.0+elipson*cos(angle))
            ft = (1.0/s1%r)*(vc**2)*log(s1%r)*(elipson*m*sin(angle))


            ! update accelerations
            s2%r2     = fr+s1%r*(s1%theta1)**2
            s2%theta2 = (ft-2.0*s1%r1*s1%theta1)/s1%r

            s2%l      = ((s1%r)**2)*(s1%theta1)
            s2%energy = 0.5*((s1%r)**2)*((s1%theta1)**2)

            s2%J      = s2%energy-omegap*s2%l
            s2%Jin    = (s2%J+1.0)/(-1.0)

        end subroutine update_acc

end module physics

