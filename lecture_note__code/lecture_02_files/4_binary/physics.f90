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
            integer :: mass,d,v
            real :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11
            real :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11
            real :: v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11
            
            mass = 10**27
            d = 10**11
            v = 10**5

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            ! a           = 3*au !TODO
            ! a3          = 10*au !TODO
            ! period1     = (((4*pi**2)/(G*(m2+m3)))*a**3)**0.5 !TODO
            ! period2     = (((4*pi**2)/(G*(m1+m2+m3)))*a3**3)**0.5 !TODO
            ! !period     = 1*365*86400 !TODO
            ! force12     = (G*m1*m2)/(8*au)**2  !TODO
            ! force13     = (G*m1*m3)/(11*au)**2  !TODO
            ! force23     = (G*m3*m2)/(3*au)**2  !TODO

            m1=msun
            m2=0.330*mass
            m3=4.87*mass	
            m4=5.97*mass	
            m5=0.073*mass	
            m6=0.642*mass	
            m7=1898*mass
            m8=568*mass
            m9=86.8*mass
            m10=102*mass
            m11=0.0146*mass

            v2=47.4*v
            v3=35.0*v
            v4=29.8*v
            v5=1.0*v
            v6=24.1*v
            v7=13.1*v
            v8=9.7*v
            v9=6.8*v
            v10=5.4*v
            v11=4.7*v

            d2=57.9*d
            d3=108.2*d
            d4=149.6*d	
            d5=0.384*d	
            d6=227.9*d	
            d7=778.6*d	
            d8=1433.5*d	
            d9=2872.5*d	
            d10=4495.1*d	
            d11=5906.4*d
            

            ! 太陽
            ! setup initial conditions of star M2
            stars(1)%mass = m1
            stars(1)%x    = 0.0 !TODO
            stars(1)%xt   = 0.0 !TODO
            stars(1)%y    = 0.0 !TODO
            stars(1)%yt   = 0.0 !TODO
            stars(1)%vx   = 0.0 !TODO
            stars(1)%vxt  = 0.0 !TODO
            stars(1)%vy   = 0.0 !TODO
            stars(1)%vyt  = 0.0 !TODO
            stars(1)%ax   = 0.0 !TODO
            stars(1)%axt  = 0.0 !TODO
            stars(1)%ay   = 0.0 !TODO
            stars(1)%ayt  = 0.0 !TODO

            ! 水星
            ! setup initial conditions of star M1
            stars(2)%mass = m2
            stars(2)%x    = d2!TODO
            stars(2)%xt   = 0.0!TODO
            stars(2)%y    = 0.0 !TODO
            stars(2)%yt   = 0.0 !TODO
            stars(2)%vx   = 0.0 !TODO
            stars(2)%vxt  = 0.0 !TODO
            stars(2)%vy   = v2 !TODO
            stars(2)%vyt  = 0.0 !TODO
            stars(2)%ax   = -((G*m1*m2)/(d2)**2)/m2 !TODO
            stars(2)%axt  = stars(2)%ax !TODO
            stars(2)%ay   = 0.0 !TODO
            stars(2)%ayt  = 0.0 !TODO

            ! 金星
            ! setup initial conditions of star M2
            stars(3)%mass = m3
            stars(3)%x    = d3 !TODO
            stars(3)%xt   = 0.0 !TODO
            stars(3)%y    = 0.0 !TODO
            stars(3)%yt   = 0.0 !TODO
            stars(3)%vx   = 0.0 !TODO
            stars(3)%vxt  = 0.0 !TODO
            stars(3)%vy   = v3 !TODO
            stars(3)%vyt  = 0.0 !TODO
            stars(3)%ax   = -((G*m1*m3)/(d3)**2)/m3 !TODO
            stars(3)%axt  = stars(3)%ax !TODO
            stars(3)%ay   = 0.0 !TODO
            stars(3)%ayt  = 0.0 !TODO

            !地球
            ! setup initial conditions of star M2
            stars(4)%mass = m4
            stars(4)%x    = d4 !TODO
            stars(4)%xt   = 0.0 !TODO
            stars(4)%y    = 0.0 !TODO
            stars(4)%yt   = 0.0 !TODO
            stars(4)%vx   = 0.0 !TODO
            stars(4)%vxt  = 0.0 !TODO
            stars(4)%vy   = v4-v5 !TODO
            stars(4)%vyt  = 0.0 !TODO
            stars(4)%ax   = -((G*m1*m4)/(d4)**2)/m4 !TODO
            stars(4)%axt  = stars(4)%ax !TODO
            stars(4)%ay   = 0.0 !TODO
            stars(4)%ayt  = 0.0 !TODO

            ! moon
            ! setup initial conditions of star M2
            stars(5)%mass = m5
            stars(5)%x    = d4+d5 !TODO
            stars(5)%xt   = 0.0 !TODO
            stars(5)%y    = 0.0 !TODO
            stars(5)%yt   = 0.0 !TODO
            stars(5)%vx   = 0.0 !TODO
            stars(5)%vxt  = 0.0 !TODO
            stars(5)%vy   = v4+v5 !TODO
            stars(5)%vyt  = 0.0 !TODO
            stars(5)%ax   = -((G*m4*m5)/(d5)**2)/m5-((G*m1*m5)/(d4+d5)**2)/m5 !TODO
            stars(5)%axt  = stars(5)%ax !TODO
            stars(5)%ay   = 0.0 !TODO
            stars(5)%ayt  = 0.0 !TODO

            ! mars
            ! setup initial conditions of star M2
            stars(6)%mass = m6
            stars(6)%x    = d6 !TODO
            stars(6)%xt   = 0.0 !TODO
            stars(6)%y    = 0.0 !TODO
            stars(6)%yt   = 0.0 !TODO
            stars(6)%vx   = 0.0 !TODO
            stars(6)%vxt  = 0.0 !TODO
            stars(6)%vy   = v6 !TODO
            stars(6)%vyt  = 0.0 !TODO
            stars(6)%ax   = -((G*m1*m6)/(d6)**2)/m6 !TODO
            stars(6)%axt  = stars(6)%ax !TODO
            stars(6)%ay   = 0.0 !TODO
            stars(6)%ayt  = 0.0 !TODO

            ! jupiter
            ! setup initial conditions of star M2
            stars(7)%mass = m7
            stars(7)%x    = d7 !TODO
            stars(7)%xt   = 0.0 !TODO
            stars(7)%y    = 0.0 !TODO
            stars(7)%yt   = 0.0 !TODO
            stars(7)%vx   = 0.0 !TODO
            stars(7)%vxt  = 0.0 !TODO
            stars(7)%vy   = v7 !TODO
            stars(7)%vyt  = 0.0 !TODO
            stars(7)%ax   = -((G*m1*m7)/(d7)**2)/m7 !TODO
            stars(7)%axt  = stars(7)%ax !TODO
            stars(7)%ay   = 0.0 !TODO
            stars(7)%ayt  = 0.0 !TODO

            ! 土星SATURN
            ! setup initial conditions of star M2
            stars(8)%mass = m8
            stars(8)%x    = d8 !TODO
            stars(8)%xt   = 0.0 !TODO
            stars(8)%y    = 0.0 !TODO
            stars(8)%yt   = 0.0 !TODO
            stars(8)%vx   = 0.0 !TODO
            stars(8)%vxt  = 0.0 !TODO
            stars(8)%vy   = v8 !TODO
            stars(8)%vyt  = 0.0 !TODO
            stars(8)%ax   = -((G*m1*m8)/(d8)**2)/m8 !TODO
            stars(8)%axt  = stars(8)%ax !TODO
            stars(8)%ay   = 0.0 !TODO
            stars(8)%ayt  = 0.0 !TODO

            ! 天文星URANUS
            ! setup initial conditions of star M2
            stars(9)%mass = m9
            stars(9)%x    = d9 !TODO
            stars(9)%xt   = 0.0 !TODO
            stars(9)%y    = 0.0 !TODO
            stars(9)%yt   = 0.0 !TODO
            stars(9)%vx   = 0.0 !TODO
            stars(9)%vxt  = 0.0 !TODO
            stars(9)%vy   = v9 !TODO
            stars(9)%vyt  = 0.0 !TODO
            stars(9)%ax   = -((G*m1*m9)/(d9)**2)/m9 !TODO
            stars(9)%axt  = stars(9)%ax !TODO
            stars(9)%ay   = 0.0 !TODO
            stars(9)%ayt  = 0.0 !TODO

                        !
            ! setup initial conditions of star M2
            stars(10)%mass = m10
            stars(10)%x    = d10 !TODO
            stars(10)%xt   = 0.0 !TODO
            stars(10)%y    = 0.0 !TODO
            stars(10)%yt   = 0.0 !TODO
            stars(10)%vx   = 0.0 !TODO
            stars(10)%vxt  = 0.0 !TODO
            stars(10)%vy   = v10 !TODO
            stars(10)%vyt  = 0.0 !TODO
            stars(10)%ax   = -((G*m1*m10)/(d10)**2)/m10 !TODO
            stars(10)%axt  = stars(10)%ax !TODO
            stars(10)%ay   = 0.0 !TODO
            stars(10)%ayt  = 0.0 !TODO

                        !
            ! setup initial conditions of star M2
            stars(11)%mass = m11
            stars(11)%x    = d11 !TODO
            stars(11)%xt   = 0.0 !TODO
            stars(11)%y    = 0.0 !TODO
            stars(11)%yt   = 0.0 !TODO
            stars(11)%vx   = 0.0 !TODO
            stars(11)%vxt  = 0.0 !TODO
            stars(11)%vy   = v11 !TODO
            stars(11)%vyt  = 0.0 !TODO
            stars(11)%ax   = -((G*m1*m11)/(d11)**2)/m11 !TODO
            stars(11)%axt  = stars(11)%ax !TODO
            stars(11)%ay   = 0.0 !TODO
            stars(11)%ayt  = 0.0 !TODO




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

                    else if (j==k) then
                        cycle
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

                ! stars(j)%axt = fxt/stars(j)%mass
                ! stars(j)%ayt = fyt/stars(j)%mass
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

