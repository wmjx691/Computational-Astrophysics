from simulation_data import N, stell_info, starinfo
from constants import G, pi, msun, mass, au, yr, sec, dis
import numpy as np
import math


def rk2(dt, yin, derive_rk2, stell):
    """
    Do one RK2 update with dt
    """
    dydt   = derive_rk2(yin, stell)
    ystell = yin+dt*dydt
    dydt   = derive_rk2(ystell, stell)
    yout   = ystell+dt*dydt
    yout   = 0.5*(yin+yout)

    return yout

def derive_rk2(yin, stell):
    """
    Input:
        yin : numpy array

    Output:
        yout : numpy array : the dydt
    """
    dydt0 = yin[3]         # dxdt  = vx
    dydt1 = yin[4]         # dydt  = vy
    dydt2 = yin[5]         # dzdt  = vz
    dydt3 = stell.ax       # dvxdt = ax
    dydt4 = stell.ay       # dvydt = ay
    dydt5 = stell.az       # dvzdt = az

    dydt  = np.array([dydt0, dydt1, dydt2, 
                      dydt3, dydt4, dydt5])
    return dydt

def euler(dt, yin1, yin2, derive, stell_acc):
    """
    Do one euler update with dt
    """
    dydt   = derive(yin2, stell_acc)
    yout   = yin1 + dt*dydt

    return yout

def derive(yin, stell_acc):
    """
    Input:
        yin : numpy array

    Output:
        yout : numpy array : the dydt
    """
    dydt0 = yin[3]         # dxdt  = vx
    dydt1 = yin[4]         # dydt  = vy
    dydt2 = yin[5]         # dzdt  = vz
    dydt3 = stell_acc[0]   # dvxdt = ax
    dydt4 = stell_acc[1]   # dvydt = ay
    dydt5 = stell_acc[2]   # dvzdt = az

    dydt  = np.array([dydt0, dydt1, dydt2, 
                      dydt3, dydt4, dydt5])
    return dydt



def update(stell, stars, step, dt):
    """
    Do one update with dt = 0.01 [sec]

    """
    yin = np.array([stell.x,
                    stell.y,
                    stell.z,
                    stell.vx,
                    stell.vy,
                    stell.vz])
    

    # #rk2
    # yout = rk2(dt, yin, derive_rk2, stell)

    # stell.x  = yout[0]
    # stell.y  = yout[1]
    # stell.z  = yout[2]
    # stell.vx = yout[3]
    # stell.vy = yout[4]
    # stell.vz = yout[5]

    # stell_acc  = update_acc(step, yin, stell.mass, stars)
    # stell.ax = stell_acc[0]
    # stell.ay = stell_acc[1]
    # stell.az = stell_acc[2]

    #rk4
    mass_stell = stell.mass

    h  =  0.5*dt
    stell_acc  = update_acc(step, yin, mass_stell, stars)
    yout2      = euler(h, yin, yin, derive, stell_acc)

    step       = step+1
    h          = 0.5*dt
    stell2_acc = update_acc(step, yout2, mass_stell, stars)
    yout3      = euler(h, yin, yout2, derive, stell2_acc)

    step       = step+1
    h  =  dt
    stell3_acc = update_acc(step, yout3, mass_stell, stars)
    yout4      = euler(h, yin, yout3, derive, stell3_acc)

    step       = step+2
    stell4_acc = update_acc(step, yout4, mass_stell, stars)

    #          x0     + dt*(vx0    + 2vx1       + 2vx2       + vx3)/6
    stell.x  = yin[0] + dt*(yin[3] + 2.*yout2[3]+ 2.*yout3[3]+ yout4[3])/6.0
    stell.y  = yin[1] + dt*(yin[4] + 2.*yout2[4]+ 2.*yout3[4]+ yout4[4])/6.0
    stell.z  = yin[2] + dt*(yin[5] + 2.*yout2[5]+ 2.*yout3[5]+ yout4[5])/6.0
    stell.vx = yin[3] + dt*(stell_acc[0] + 2.*stell2_acc[0]+ 2.*stell3_acc[0]+ stell4_acc[0])/6.0
    stell.vy = yin[4] + dt*(stell_acc[1] + 2.*stell2_acc[1]+ 2.*stell3_acc[1]+ stell4_acc[1])/6.0
    stell.vz = yin[5] + dt*(stell_acc[2] + 2.*stell2_acc[2]+ 2.*stell3_acc[2]+ stell4_acc[2])/6.0


    asave = np.array([stell.x,
                    stell.y,
                    stell.z,
                    stell.vx,
                    stell.vy,
                    stell.vz])

    stell_save_acc = update_acc(step, asave, mass_stell, stars)

    stell.ax    = stell_save_acc[0]
    stell.ay    = stell_save_acc[1]
    stell.az    = stell_save_acc[2]

    return stell

# update the stell new ax,ay,az    
def update_acc(step, yin, mass_stell, stars):
    fx=0
    fy=0
    fz=0

    x     = yin[0] - stars[9].x[step]
    y     = yin[1] - stars[9].y[step]
    z     = yin[2] - stars[9].z[step]

    rsq   = x**2 + y**2 + z**2
    rflat = x**2 + y**2
    angle = math.atan2(y,x)
    phi   = math.atan2(z,rflat**0.5)
    force = G * mass_stell * stars[9].mass /rsq


    fx = fx - force*math.cos(phi)*math.cos(angle)
    fy = fy - force*math.cos(phi)*math.sin(angle)
    fz = fz - force*math.sin(phi)

    stell_acc = np.array([  fx/mass_stell,
                            fy/mass_stell,
                            fz/mass_stell])

    return stell_acc

# getting the angle relation of this vector
def get_angle(x ,y, z):

    rflat = x**2 + y**2
    angle = math.atan2(y,x)
    phi   = math.atan2(z,rflat**0.5)

    vector_transform  =  np.array([
                                    math.cos(phi)*math.cos(angle),
                                    math.cos(phi)*math.sin(angle),
                                    math.sin(phi)
    ])

    return vector_transform


# getting escape vel
def soi(escape_mass, distance):
    escape_vel = (2.*G*escape_mass/distance)**0.5

    return escape_vel

# if __name__=='__main__':

    mass_star   = np.array([  0.330 *mass,   # 0Mercury
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

    radius_star = np.array([  4879,
                            12104,
                            12756,
                            6792,
                            142984,
                            120536,
                            51118,
                            49528,
                            2370,
                            1392700
    ])
    radius_star = radius_star/2

    Escape_Vel  = np.array([4.3, 10.4, 11.2, 5.0, 59.5, 35.5, 21.3, 23.5, 1.3, 615.])
    Escape_Vel  = Escape_Vel*sec

    dis_stars   = np.array([57.9,108.2,149.6,227.9,778.6,1433.5,2872.5,4495.1,5906.4,0.])
    dis_stars   = dis_stars*dis

    # calculate sphere of influence of each planet(soi)


#     m1 = 1.0 * msun
#     m2 = 2.0 * msun
#     m3 = 1.0 * msun
#     a           = 3*au
#     a3          = 10*au
#     period1     = (((4*pi**2)/(G*(m1+m2)))*a**3)**0.5
#     period2     = (((4*pi**2)/(G*(m1+m2+m3)))*a3**3)**0.5
#     force12     = (G*m1*m2)/(3*au)**2
#     force13     = (G*m1*m3)/(8*au)**2
#     force23     = (G*m3*m2)/(11*au)**2

    
#     star1=np.array([1,        #id
#                     m1,       #mass
#                     -2.0*au,  #x
#                     0,        #y
#                     0,        #vx
#                     (2.0*pi*m2/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2),
#                     (-force13+force12)/m1, 
#                     0         #ay
#     ])

#     star2=np.array([2,        #id
#                     m2,       #mass
#                     au,       #x
#                     0,        #y
#                     0,        #vx
#                     (-2.0*pi*m1/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2),
#                     -(force12+force23)/m2, 
#                     0         #ay
#     ])

#     star3=np.array([3,        #id
#                     m3,       #mass
#                     -10.0*au,  #x
#                     0,        #y
#                     0,        #vx
#                     (2.0*pi*(m1+m2)/(m1+m2+m3)*a3/period2),
#                     (force23+force13)/m3, 
#                     0         #ay
#     ])

#     starnp=np.array([star1, star2, star3])

#     for n in np.arange(100):
#         update(starnp, dt=0.01)
#         # print(point.posx,point.posy)