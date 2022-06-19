#binary
from simulation_data import N, stell_info
from constants import G, pi, msun, mass, au, yr
from physics import update, update_acc
from output import record
import numpy as np

from inputfile import ephem, set_starinfo


# initial setup
dt       = 1
t_launch = "1977-08-20 00:00"   #lanuch from earth
t_arrive = "1988-08-03 00:00"   #arrive to   Neptune

# 2dt is prepare to do rk4
position  , velcoity  , time_duration   = ephem(t_launch, t_arrive, 0.5*dt)


# output intervals to record trajectory
step     = 0             #start from 0th step
noutput  = 1000
interval = (len(time_duration)/dt)/noutput


# initialize the stars
mass_star = np.array([  0.330 *mass,   # 0Mercury
                        4.87  *mass,   # 1Venus
                        5.97  *mass,   # 2Earth
                        0.642 *mass,   # 3Mars
                        1898. *mass,   # 4Jupiter
                        568.  *mass,   # 5Saturn
                        86.8  *mass,   # 6Uranus
                        102.  *mass,   # 7Neptune
                        0.0146*mass,   # 8Pluto
                               msun,   # 9sun
])


# calculate sphere of influence of each planet



#stars array has stars[i]:(idindex, mass, x, y)
stars   = set_starinfo(mass_star, position  , velcoity)


# #create the record of each star
# position_o, velcoity_o, time_duration_o = ephem(t_launch, t_arrive, dt)
# stars_o = set_starinfo(mass_star, position_o, velcoity_o)
# #record each planet in solar system
# record(stars_o, time_duration_o, first, planet)

#output each stars info
first=True
planet=True


# initialize the stellite
mass_stell = 721.9  #kg

stell_setup = np.array([stars[2].x[0],   #x
                        stars[2].y[0],   #y
                        stars[2].z[0],   #z
                        stars[2].vx[0],  #vx
                        stars[2].vy[0],  #vy
                        stars[2].vz[0]  #vz
])

stell_acc = update_acc(step, stell_setup, mass_stell, stars)

stell = stell_info(         10,  #idindex
                    mass_stell,  #mass
                stell_setup[0],   #x
                stell_setup[1],   #y
                stell_setup[2],   #z
                stell_setup[3],  #vx
                stell_setup[4],  #vy
                stell_setup[5],  #vz
                stell_acc[0],    #ax
                stell_acc[1],    #ay
                stell_acc[2]     #az
)

# # # test for the earth 2d model
# r        = 5.2*au
# period   = (((4.0*pi**2)/(G*(msun+stars[4].mass)))*r**3)**0.5
# force    = (G*stars[4].mass*msun)/(r)**2

# stell = stell_info(         10,  #idindex
#                     stars[4].mass,  #mass
#                 r,   #x
#                 0,   #y
#                 0,   #z
#                 0,   #vx
#                 (2.0*pi*r/period),  #vy
#                 0,  #vz
#                 -force/stars[4].mass,    #ax
#                 0,    #ay
#                 0     #az
# )



planet = False

for time in time_duration:
    if (step == 7998/dt):
        break

    if ((time-2443376.5)%1. != 0):
        step = step+1
        continue
    
    stell = update(stell, stars, step, dt)

    first = record(stell, time, first, planet)

    step = step+1

print("Done for simulation!")