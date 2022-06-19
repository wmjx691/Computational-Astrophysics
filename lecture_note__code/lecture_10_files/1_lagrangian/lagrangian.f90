      program lh1
!
!     This Lagrangian code follows the adiabatic expansion 
!     of a hot spherical bubble in a uniform ambient medium.
!     The present setup is ajusted to spherical geometry.
!     The results may be plotted in Supermomgo environment with the help
!     of commands contained in files pltmod and plte.
!     The code can also work in Cartesian geometry. To switch from 
!     spherical to Cartesian geometry go to subroutine inicond and
!     set cartesian = .true.
!
      use IO
      use Simulation_data
      implicit none

      integer :: ix, nsteps, iout
      real :: ekin, emax, epsi, ethe, etot, etot0
      real :: pmax, rhomax, umax, wmax, dmamb,dmej, rhoamb, rhoej
!
!                                  open the input-data file
!
      open(15,file='input.dat',form='formatted')
!
!                                  open the file in which 
!                                  the final model will be stored 
!
      open(16,file='lh1.out',form='formatted')
!
!                                  open the file in which values
!                                  of thermal, kinetic and total
!                                  energy will be stored at every
!                                  10th time-step
!
      open(17,file='lh2.out',form='formatted')
!
      read(15,*)
      read(15,*)nsteps
      read(15,*)tmax   
      read(15,*)efac   
      read(15,*)il
      read(15,*)q
      read(15,*)fmratio
!
! set initial conditions
!
      call initial()

!======================================================
!
!
! begin loop over time steps
!
!
!======================================================

      iout  = 0
      istep = 1

      do while((t .le. tmax) .and. (istep .le. nsteps))

!
! determine time step
!
        call find_dt()
!
! update velocities
!
        do ix=2,nx
          u(ix) = u(ix)-a(ix)*(p(ix)-p(ix-1))*dt/dm(ix) &
                -0.5*(w(ix  )*(3.*ak12(ix  )-a(ix)) &
                     -w(ix-1)*(3.*ak12(ix-1)-a(ix))) &
                    *dt/dm(ix)
        end do
        if (cartesian) then
          u(1) = u(2)
        else
          u(1) = 0.
        end if
!
! update radii, surfaces and volumes
!
        do ix=1,nx
          rold(ix) = r(ix)
        end do
        do ix=1,nx
          r(ix)    = rold(ix)+u(ix)*dt12
        end do
        do ix=1,nx-1
          dr12(ix) = r(ix+1)-r(ix)
        end do
        if (cartesian) then
          do ix=1,nx
            v(ix) = r(ix)
          end do
        else
          do ix=1,nx
            at12(ix) = 4.*pi*(0.5*(r(ix)+rold(ix)))**2
            a   (ix) = 4.*pi*r(ix)**2 
            v   (ix) = 4./3.*pi*r(ix)**3 
          end do
          do ix=1,nx-1
            ak12(ix) = 0.5*(at12(ix+1)+at12(ix))
          end do
        end if
!
! update densities
!
        do ix=1,nx-1
          rho(ix)  = dm12(ix)/(v(ix+1)-v(ix))
        end do
          rho(nx)  = rho(nx-1)
!
! artificial viscosity
!
        do ix=1,nx-1
          w(ix)  =-q**2*rho(ix)*abs(u(ix+1)-u(ix)) &
                   *(u(ix+1)*(1.-at12(ix+1)/3./ak12(ix)) &
                    -u(ix  )*(1.-at12(ix  )/3./ak12(ix)))
        end do
        do ix=1,nx-1
          if (u(ix+1).gt.u(ix)) w(ix)=0.
        end do
!
! update internal energies and pressures
!
        do ix=1,nx-1
          aux(ix)  = eps(ix)-p(ix) &
                   *( at12(ix+1)*u(ix+1) &
                     -at12(ix  )*u(ix  ))*dt12/dm12(ix)
        end do
        do ix=1,nx-1
          p(ix)    = 0.5*(p(ix)+(gamma-1.)*rho(ix)*aux(ix))
        end do
        do ix=1,nx-1
          eps(ix)  = eps(ix)-p(ix) &
                   *( at12(ix+1)*u(ix+1) &
                     -at12(ix  )*u(ix  ))*dt12/dm12(ix)
        end do

!
! contribution from artificial viscosity
!

        do ix=1,nx-1
          eps(ix)  = eps(ix)-0.5*w(ix)*dt12/dm12(ix) &
                   *(u(ix+1)*(3.*ak12(ix)-at12(ix+1)) &
                    -u(ix  )*(3.*ak12(ix)-at12(ix  )))
        end do 
!
        do ix=1,nx-1
          p(ix)    = (gamma-1.)*rho(ix)*eps(ix)
        end do
        p  (nx) = p  (nx-1)
        eps(nx) = eps(nx-1)
!
!                                 end time step
!
        t = t+dt12  
!
!                                 check energy conservstion
!
        ethe = 0.
        ekin = 0.
        do ix=2,nx-1
          ethe = ethe + (eps(ix))*dm(ix)
          ekin = ekin + 0.5*(0.5*(u(ix+1)+u(ix)))**2*dm(ix)
        end do
        etot = ethe + ekin
        if (istep.eq.1) etot0 = etot
        etot = etot/etot0
        ethe = ethe/etot0
        ekin = ekin/etot0
        if (mod(istep,200).eq.0) then
            write(17,200)istep,t,etot, ethe, ekin
            call output(iout,t)
            iout = iout + 1
        endif
        if (mod(istep,100).eq.0) write(* ,100)istep,t,etot, ethe, ekin
!
!
!                                      end loop over time steps
!
        istep = istep + 1
      end do
!
!
!                                      final printout
!

      umax=0.
      rhomax=0.
      pmax=0.
      emax=0.
      wmax=0.

      do ix=1,nx
         if (u(ix).gt.0) umax  = max(u(ix),umax)
         rhomax = max(rho(ix),rhomax)
         pmax   = max(p(ix),pmax)
         emax   = max(eps(ix),emax)
         wmax   = max(w(ix),wmax)
      end do

      epsi=1.e-20

      do ix=1,nx
         u(ix)   = u(ix)/(epsi+umax)
         rho(ix) = rho(ix)/(epsi+rhomax)
         p(ix)   = p(ix)/(epsi+pmax)
         eps(ix) = eps(ix)/(epsi+emax)
         w(ix)   = w(ix)/(epsi+wmax)
      end do

      write(16,101)(ix,fm(ix),r(ix)/r(nx),u(ix), &
            rho(ix),p(ix),eps(ix),w(ix),ix=1,nx)
  100 format(1x,i5,'; t:',1pe10.2,';  etot:',1pe10.2, &
                ';  eth:',1pe10.2,';  ekin:',1pe10.2)
  101 format(1x,i5,1p7e12.4) 
  200 format(1x,i5,1p4e11.3)
!
      stop
      end program
!=================================================================
!
!
!
!
!=================================================================
      subroutine initial
      use Simulation_data
      implicit none

      integer :: ix
      real :: rhoej, dmamb, dmej, rhoamb
!
      cartesian = .false.
!
!     Courant factor, must be smaller than 1
      cflfactor = 0.1
!
!     diffusion factor, must be greater that 1
      dfactor   = 1.5
      dfactor   = dfactor**2
!
      gamma     = 5.0/3.0
!
      pi        = asin(1.0)*2.
!
!                                         high pressure ejecta
!
      rhoej  =  ((float(nx)/float(il))**3-1.)*fmratio
      dmej   =  4./3.*pi*(float(il)/float(nx))**3*rhoej/float(il)
      do ix = 1,il
        dm12(ix) = dmej
        rho (ix) = rhoej
        eps (ix) = efac
        p   (ix) = (gamma-1.)*rho(ix)*eps(ix)
        u   (ix) = 0.
      end do
!     
!                                         low pressure ambient medium
!
      rhoamb = 1.
      dmamb  = 4./3.*pi*(1.-(float(il)/float(nx))**3)/float(nx-il)*rhoamb
      do ix=il+1,nx
        dm12(ix) = dmamb
        rho (ix) = rhoamb
        eps (ix) = 1.
        p   (ix) = (gamma-1.)*rho(ix)*eps(ix)
        u   (ix) = 0. 
      end do
!
      fm(1) = dm12(1) 
      do ix=2,nx
        fm (ix) = fm(ix-1)+dm12(ix)
      end do
      do ix=2,nx
        dm(ix) = 0.5*(dm12(ix)+dm12(ix-1))
      end do
!
!
      r   (1) = 0.
      v   (1) = 0.
      do ix=2,nx
        v(ix) = v(ix-1)+dm12(ix-1)/rho(ix-1)
        if (cartesian) then
          r(ix)    = v(ix)
          a(ix)    = 1.
          at12(ix) = 1.
        else
          r(ix) = (v(ix)/(4./3.*pi))**(1./3.)
          a(ix) = 4.*pi*r(ix)**2
        end if
      end do
      do ix=1,nx-1
        dr12(ix) = r(ix+1)-r(ix)
      end do
!
      t    = 0.
!
      return
      end
!=================================================================
!
!
!
!=================================================================
      subroutine find_dt
      use Simulation_data
      implicit none
      real :: dtc, dtd
      integer :: ix
!
      dtc = 1.e30
      do ix=1,nx-1
        dtc = min(dtc,dr12(ix) &
               /(abs(u(ix))+sqrt(gamma*eps(ix))))
      end do
!
      dtc  = cflfactor*dtc
      if (t+dtc.gt.tmax) dtc=tmax-t 
!
!     diffusion limit
!
      dtd = 1.e-30
      do ix=1,nx-1
        dtd = max(dtd,abs(at12(ix+1)*u(ix+1)-at12(ix)*u(ix)) &
                       /(   v(ix+1)        -v(ix)         ))
      end do
      dtd = 0.5/dtd/dfactor
!
      dtc = min(dtc,dtd)
!
      dt   = 0.5*(dt12+dtc)
      dt12 = dtc
!
      return
      end
