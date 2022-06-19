from simulation_data import N
from constants import *
import numpy as np
import math

def rk2(dt, yin, derive, starnp, i):
    """
    Do one RK2 update with dt
    """
    dydt  = derive(yin, starnp, i)
    ystar = yin+dt*dydt
    dydt  = derive(ystar, starnp, i)
    yout  = ystar+dt*dydt
    yout  = 0.5*(yin+yout)
    return yout

def derive(yin, starnp, i):
    """
    Input:
        yin : numpy array

    Output:
        yout : numpy array : the dydt
    """
    dydt0 = yin[2]      # dxdt = vx
    dydt1 = yin[3]      # dydt = vy
    dydt2 = starnp[i][6]      # dvxdt = ax
    dydt3 = starnp[i][7]      # dvydt = ay

    dydt = np.array([dydt0, dydt1, 
                    dydt2, dydt3])
    return dydt



def update(starnp, dt):
    """
    Do one update with dt = 0.01 [sec]

    """
    for i in np.arange(0, N):
        yin = np.array([starnp[i][2],
                        starnp[i][3],
                        starnp[i][4],
                        starnp[i][5]])

        yout = rk2(dt, yin, derive, starnp, i)

        starnp[i][2]=yout[0]
        starnp[i][3]=yout[1]
        starnp[i][4]=yout[2]
        starnp[i][5]=yout[3]

    
    for j in np.arange(0, N):
        fx=0
        fy=0
    
        for k in np.arange(0, N):
            if j==k:
                continue
            if j<k:
                x =starnp[j][2]-starnp[k][2]
                y =starnp[j][3]-starnp[k][3]      
            elif j>k:
                x =-(starnp[j][2]-starnp[k][2])
                y =-(starnp[j][3]-starnp[k][3])

            rsq = x**2 + y**2
            angle = math.atan2(y,x)
            force = G *starnp[j][1] *starnp[k][1] /rsq

            if j<k:
                fx=fx-force*math.cos(angle)
                fy=fy-force*math.sin(angle)
            elif j>k:
                fx=fx+force*math.cos(angle)
                fy=fy+force*math.sin(angle)
                
        starnp[j][6] = fx/starnp[j][1]
        starnp[j][7] = fy/starnp[j][1]


    return starnp

        
def update_acc(dt,starnp):

    return


# if __name__=='__main__':
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