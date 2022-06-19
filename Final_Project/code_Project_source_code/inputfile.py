from constants import G, pi, msun, mass, au, yr
from simulation_data import starinfo, N

from astropy.time import Time
import numpy as np
from jplephem.spk import SPK



def ephem(t_launch, t_arrive, dt):

    t_launch = Time(t_launch)
    t_arrive = Time(t_arrive)
    
    
    time_duration = np.arange(t_launch.jd,t_arrive.jd,dt)
    time_test = np.arange(t_launch.jd,t_launch.jd+365,dt)
    
    # call jpl
    kernel = SPK.open('de421.bsp')
    
    
    ####change back to time_duration
    position = np.zeros((N,3,len(time_duration)))
    velocity = np.zeros((N,3,len(time_duration)))
    
    for i in np.arange(0, N):
        position[i][:][:],velocity[i][:][:] = kernel[0,i+1].compute_and_differentiate(time_duration)
                
    
    # position[which stars?][x,y,z?][time?]
    # velocity[which stars?][x,y,z?][time?]
    return position, velocity, time_duration



def set_starinfo(mass_star, position, velcoity, soi):
    stars   = np.zeros(N, dtype=np.dtype(starinfo))
    
    for i in np.arange(0, N):
        #star   = starinfo(idindex,         mass,
        #                                      x,                 y,                 z)
        stars[i]= starinfo(      i, mass_star[i],
                               position[i][0][:], position[i][1][:], position[i][2][:],
                               velcoity[i][0][:], velcoity[i][1][:], velcoity[i][2][:],
                               soi[i])

        # 0Mercury
        # 1Venus
        # 2Earth
        # 3Mars
        # 4Jupiter
        # 5Saturn
        # 6Uranus
        # 7Neptune
        # 8Pluto
        # 9sun 

    return  stars
 
   
#if __name__=='__main__':

    dt       = 1    
    t_launch = "1977-08-20 12:00"   #lanuch from earth
    t_arrive = "1987-08-20 12:00"   #arrive to   Neptune
    
    position, velocity, time_duration = ephem(t_launch, t_arrive, dt)
    
    mass_star = np.array([   msun,   # 0sun
                      0.330 *mass,   # 1Mercury
                      4.87  *mass,   # 2Venus
                      5.97  *mass,	 # 3Earth
                      0.642 *mass,	 # 4Mars
                      1898. *mass,   # 5Jupiter
                      568.  *mass,   # 6Saturn
                      86.8  *mass,   # 7Uranus
                      102.  *mass,   # 8Neptune
                      0.0146*mass,   # 9Pluto
    ])
    
    stars = set_starinfo(mass_star, position, velocity)
    
    print(stars[2].get_loc(0))
