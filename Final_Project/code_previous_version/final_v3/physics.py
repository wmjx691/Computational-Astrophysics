from simulation_data import N, stell_info, speedupinfomation
from constants import G, sec, fuel_acc
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



def update(stell, stars, step, dt, speedupinfo):
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

    h  =  0.5*dt
    stell_acc  = update_acc(step, yin, stell, stars, speedupinfo)
    yout2      = euler(h, yin, yin, derive, stell_acc)

    step       = step+1
    h          = 0.5*dt
    stell2_acc = update_acc(step, yout2, stell, stars, speedupinfo)
    yout3      = euler(h, yin, yout2, derive, stell2_acc)

    step       = step+1
    h  =  dt
    stell3_acc = update_acc(step, yout3, stell, stars, speedupinfo)
    yout4      = euler(h, yin, yout3, derive, stell3_acc)

    step       = step+2
    stell4_acc = update_acc(step, yout4, stell, stars, speedupinfo)

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

    stell_save_acc = update_acc(step, asave, stell, stars, speedupinfo)

    stell.ax    = stell_save_acc[0]
    stell.ay    = stell_save_acc[1]
    stell.az    = stell_save_acc[2]

    return stell

# update the stell new ax,ay,az    
def update_acc(step, yin, stell, stars, speedupinfo):
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
    force = G * stell.mass * stars[9].mass /rsq

    fx = - force*math.cos(phi)*math.cos(angle)
    fy = - force*math.cos(phi)*math.sin(angle)
    fz = - force*math.sin(phi)

    # for planet
    # if (soi_result.soiornot == True):
    xp     = yin[0] - stars[4].x[step]
    yp     = yin[1] - stars[4].y[step]
    zp     = yin[2] - stars[4].z[step]

    rsqp   = xp**2 + yp**2 + zp**2
    rflatp = xp**2 + yp**2
    anglep = math.atan2(yp,xp)
    phip   = math.atan2(zp,rflatp**0.5)
    forcep = G * stell.mass * stars[4].mass /rsqp

    fx = fx - forcep*math.cos(phip)*math.cos(anglep)
    fy = fy - forcep*math.cos(phip)*math.sin(anglep)
    fz = fz - forcep*math.sin(phip)

    stell_acc = np.array([  fx/stell.mass,
                            fy/stell.mass,
                            fz/stell.mass])

    if (speedupinfo.speed_ox):
        speedup_acc_vector = fuel_acc * speedupinfo.direction
        speedup_acc_vector = speedup_acc_vector * speedupinfo.dt_ratio

        stell_acc[0] = stell_acc[0] + speedup_acc_vector[0]
        stell_acc[1] = stell_acc[1] + speedup_acc_vector[1]
        stell_acc[2] = stell_acc[2] + speedup_acc_vector[2]

    return stell_acc

# getting the angle relation of this vector
def get_angle(x, y, z):

    rflat = x**2 + y**2
    angle = math.atan2(y,x)
    phi   = math.atan2(z,rflat**0.5)

    vector_transform  =  np.array([
                                    math.cos(phi)*math.cos(angle),
                                    math.cos(phi)*math.sin(angle),
                                    math.sin(phi)
    ])

    return vector_transform

# distance between stell and planet
def distance(stell, stars, idindex, step):
    x   = stell.x - stars[idindex].x[step]
    y   = stell.y - stars[idindex].y[step]
    z   = stell.z - stars[idindex].z[step]
    dis = (x**2 + y**2 + z**2)**0.5

    return dis
    

# determine that whehter "enter the soi" of this planet
def speedup_deter(stell, stars, speedupinfo, step):
    dis = distance(stell, stars, speedupinfo.idindex, step)

    if (dis <= stars[speedupinfo.idindex].rsoi) :
        speedupinfo.soi_ox = True

    else:
        speedupinfo.soi_ox = False

    return speedupinfo


# determine closest or not
def closest(stell, stars, speedupinfo, step, passrsq):
    dis = distance(stell, stars, speedupinfo.idindex, step)

    if (dis > passrsq) :
        closest_result = True
        passrsq        = passrsq

    else:
        closest_result = False
        passrsq  = dis


    return closest_result, passrsq



# getting escape vel
def escape_vel(escape_mass, distance):
    escape_vel = (2.*G*escape_mass/distance)**0.5

    return escape_vel

#if __name__=='__main__':

    fx=0

    print(fx)

    if (True):
        fx = fx - 1

    print(fx)