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
            m2 = 2.0 * msun

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            separation = 0 !TODO
            period     = 0 !TODO
            force      = 0 !TODO

            !
            ! setup initial conditions of star M1
            stars(1)%mass = m1
            stars(1)%x    = 0.0 !TODO
            stars(1)%y    = 0.0 !TODO
            stars(1)%vx   = 0.0 !TODO
            stars(1)%vy   = 0.0 !TODO
            stars(1)%ax   = 0.0 !TODO
            stars(1)%ay   = 0.0 !TODO

            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = 0.0 !TODO
            stars(2)%y    = 0.0 !TODO
            stars(2)%vx   = 0.0 !TODO
            stars(2)%vy   = 0.0 !TODO
            stars(2)%ax   = 0.0 !TODO
            stars(1)%ay   = 0.0 !TODO


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

            ! update position to t = t + dt
            !TODO

            ! update velocity to t = t + dt
            !TODO

            ! update accelerations to t = t + dt
            !TODO

            return
        end subroutine update

end module physics

