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
            real :: a, a3, m1, m2, m3, force12, force13, force23, period1, period2

            m1 = 1.0 * msun
            m2 = 2.0 * msun
            m3 = 1.0 * msun

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            a           = 3*au !TODO
            a3          = 10*au !TODO
            period1     = (((4*pi**2)/(G*(m2+m1)))*a**3)**0.5 !TODO
            period2     = (((4*pi**2)/(G*(m1+m2+m3)))*a3**3)**0.5 !TODO
            !period     = 1*365*86400 !TODO
            force12     = (G*m1*m2)/(3*au)**2  !TODO
            force13     = (G*m1*m3)/(8*au)**2  !TODO
            force23     = (G*m3*m2)/(11*au)**2  !TODO

            !
            ! setup initial conditions of star M2
            stars(3)%mass = m3
            stars(3)%x    = -10.0*au !TODO
            stars(3)%xt   = 0.0!TODO
            stars(3)%y    = 0.0 !TODO
            stars(3)%yt   = 0.0 !TODO
            stars(3)%vx   = 0.0 !TODO
            stars(3)%vxt  = 0.0 !TODO
            stars(3)%vy   = (2.0*pi*(m1+m2)/(m1+m2+m3)*a3/period2) !TODO
            stars(3)%vyt  = 0.0 !TODO
            stars(3)%ax   = (force23+force13)/m3 !TODO
            !stars(3)%axt  = stars(3)%ax !TODO
            stars(3)%ay   = 0.0 !TODO
            !stars(3)%ayt  = 0.0 !TODO

            !
            ! setup initial conditions of star M1
            stars(2)%mass = m2
            stars(2)%x    = au!TODO
            stars(2)%xt   = 0.0!TODO
            stars(2)%y    = 0.0 !TODO
            stars(2)%yt   = 0.0 !TODO
            stars(2)%vx   = 0.0 !TODO
            stars(2)%vxt  = 0.0 !TODO
            stars(2)%vy   = -(2.0*pi*m3/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2) !TODO
            stars(2)%vyt  = 0.0 !TODO
            stars(2)%ax   = -(force12+force23)/m2 !TODO
            !stars(2)%axt  = stars(2)%ax !TODO
            stars(2)%ay   = 0.0 !TODO
            !stars(2)%ayt  = 0.0 !TODO

            !
            ! setup initial conditions of star M2
            stars(1)%mass = m1
            stars(1)%x    = -2*au !TODO
            stars(1)%xt   = 0.0 !TODO
            stars(1)%y    = 0.0 !TODO
            stars(1)%yt   = 0.0 !TODO
            stars(1)%vx   = 0.0 !TODO
            stars(1)%vxt  = 0.0 !TODO
            stars(1)%vy   = (2.0*pi*m2/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2) !TODO
            stars(1)%vyt  = 0.0 !TODO
            stars(1)%ax   = (force12-force13)/m1 !TODO
            !stars(1)%axt  = stars(1)%ax !TODO
            stars(1)%ay   = 0.0 !TODO
            !stars(1)%ayt  = 0.0 !TODO


        end subroutine initial

        subroutine update(time,dt)
            use constants
            implicit none
            real, intent(in)  :: time, dt
            integer :: i,j,k
            ! real    :: x12, x13, x23, y12, y13, y23, rsq12, rsq13, rsq23
            ! real    :: fx1, fx2, fx3, fy1, fy2 ,fy3
            ! real    :: radius, force12, force13, force23, angle12, angle13, angle23

            real    :: x, y, rsq
            real    :: fx, fy
            real    :: radius, force, angle
            real    :: xt, yt, rsqt
            real    :: fxt, fyt
            real    :: radiust, forcet, anglet
            !
            ! In this example we use a first order scheme (Euler method)
            ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
            ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
            !
            ! the same approximation can be applied to dv/dt = a
            !

            do i =1, N
                ! update pos/vel/acc/time
                ! 先更新 位置pos, 在更新 速度vel
                stars(i)%xt = stars(i)%x +stars(i)%vx *dt
                stars(i)%x = 0.5*(stars(i)%x +(stars(i)%xt +stars(i)%vxt *dt))
                stars(i)%yt = stars(i)%y +stars(i)%vy *dt
                stars(i)%y = 0.5*(stars(i)%y +(stars(i)%yt +stars(i)%vyt *dt))
                stars(i)%vxt = stars(i)%vx +stars(i)%ax *dt
                stars(i)%vx = 0.5*(stars(i)%vx +(stars(i)%vxt +stars(i)%ax *dt))
                stars(i)%vyt = stars(i)%vy +stars(i)%ay *dt
                stars(i)%vy = 0.5*(stars(i)%vy +(stars(i)%vyt +stars(i)%ay *dt))
            enddo


            do j=1, N
                fx=0
                fy=0
                fxt=0
                fyt=0
                do k=1 ,N
                    if (j<k) then
                        x = stars(j)%x - stars(k)%x
                        y = stars(j)%y - stars(k)%y
                    else if (j>k) then
                        x = -(stars(j)%x - stars(k)%x)
                        y = -(stars(j)%y - stars(k)%y)
                    endif

                    rsq = x*x + y*y
                    angle = atan(y,x)
                    force = G *stars(j)%mass *stars(k)%mass /rsq

                    if (j<k) then
                        fx=fx-force*cos(angle)
                        fy=fy-force*sin(angle)
                    else if (j>k) then
                        fx=fx+force*cos(angle)
                        fy=fy+force*sin(angle)
                    endif

                enddo
                stars(j)%ax = fx/stars(j)%mass
                stars(j)%ay = fy/stars(j)%mass

                !stars(j)%axt = fxt/stars(j)%mass
                !stars(j)%ayt = fyt/stars(j)%mass
            enddo

            ! ! m1 m2 之間的距離
            ! x12 = stars(1)%x - stars(2)%x
            ! y12 = stars(1)%y - stars(2)%y
            ! x13 = stars(3)%x - stars(1)%x
            ! y13 = stars(3)%y - stars(1)%y
            ! x23 = stars(3)%x - stars(2)%x
            ! y23 = stars(3)%y - stars(2)%y

            ! rsq12 = x12*x12 + y12*y12
            ! rsq13 = x13*x13 + y13*y13
            ! rsq23 = x23*x23 + y23*y23

            ! angle12 = atan2(y12,x12)
            ! angle13 = atan2(y13,x13)
            ! angle23 = atan2(y23,x23)

            ! force12 = G *stars(1)%mass *stars(2)%mass /rsq12
            ! force13 = G *stars(1)%mass *stars(3)%mass /rsq13
            ! force23 = G *stars(2)%mass *stars(3)%mass /rsq23

            ! fx1=-force12*cos(angle12)+force13*cos(angle13)
            ! fy1=-force12*sin(angle12)+force13*sin(angle13)
            ! fx2=force12*cos(angle12)+force23*cos(angle23)
            ! fy2=force12*sin(angle12)+force23*sin(angle23)
            ! fx3=-force13*cos(angle13)-force23*cos(angle23)
            ! fy3=-force13*sin(angle13)-force23*sin(angle23)

            ! stars(1)%ax = fx1/stars(1)%mass
            ! stars(1)%ay = fy1/stars(1)%mass

            ! stars(2)%ax = fx2/stars(2)%mass
            ! stars(2)%ay = fy2/stars(2)%mass

            ! stars(3)%ax = fx3/stars(3)%mass
            ! stars(3)%ay = fy3/stars(3)%mass
            ! update position to t = t + dt
            !TODO

            ! update velocity to t = t + dt
            !TODO

            ! update accelerations to t = t + dt
            !TODO

            return
        end subroutine update

end module physics

