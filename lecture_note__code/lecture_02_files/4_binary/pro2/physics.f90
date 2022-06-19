!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        subroutine initial()
            !
            ! setup initial conditions of each stars
            ! in this example we only have two stars
            !
            use constants, only : au, msun, pi, G
            implicit none
            integer :: i
            real :: m1, m2, force

            m1 = 1.0 * msun
            m2 = 5.97219e27

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            separation = 1*au !TODO
            period     = (((4*pi**2)/(G*(m1+m2)))*separation**3)**0.5 !TODO
            !period     = 1*365*86400 !TODO
            force      = (G*m1*m2)/(separation)**2 !TODO

            !
            ! setup initial conditions of star M1
            stars(1)%mass = m1
            stars(1)%x    = -(m2/(m1+m2))*separation!TODO
            stars(1)%y    = 0.0 !TODO
            stars(1)%vx   = 0.0 !TODO
            stars(1)%vy   = 2.0*pi*m2/(m1+m2)*separation/period !TODO
            stars(1)%ax   = force/m1 !TODO
            stars(1)%ay   = 0.0 !TODO

            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = (m1/(m1+m2))*separation !TODO
            stars(2)%y    = 0.0 !TODO
            stars(2)%vx   = 0.0 !TODO
            stars(2)%vy   = 1.25*(-2.0*pi*m1/(m1+m2)*separation/period) !TODO
            stars(2)%ax   = -force/m2 !TODO
            stars(2)%ay   = 0.0 !TODO


        end subroutine initial

        subroutine update(time,dt)
            use constants
            implicit none
            real, intent(in)  :: time, dt
            integer :: i,j
            real    :: x, y, rsq, fx, fy
            real    :: radius, force, angle

            !
            ! In this example we use a first order scheme (Euler method)
            ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
            ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
            !
            ! the same approximation can be applied to dv/dt = a
            !

            do i =1, 2
                ! update pos/vel/acc/time
                ! 先更新 位置pos, 在更新 速度vel
                stars(i)%x = stars(i)%x +stars(i)%vx *dt
                stars(i)%y = stars(i)%y +stars(i)%vy *dt
                stars(i)%vx = stars(i)%vx +stars(i)%ax *dt
                stars(i)%vy = stars(i)%vy +stars(i)%ay *dt
            enddo

            ! m1 m2 之間的距離
            x = stars(1)%x - stars(2)%x
            y = stars(1)%y - stars(2)%y
            rsq = x*x + y*y
            angle = atan2(y,x)
            force = G *stars(1)%mass *stars(2)%mass /rsq
            fx=force*cos(angle)
            fy=force*sin(angle)

            stars(1)%ax = -fx/stars(1)%mass
            stars(1)%ay = -fy/stars(1)%mass

            stars(2)%ax = fx/stars(2)%mass
            stars(2)%ay = fy/stars(2)%mass

            ! update position to t = t + dt
            !TODO

            ! update velocity to t = t + dt
            !TODO

            ! update accelerations to t = t + dt
            !TODO

            return
        end subroutine update

end module physics

