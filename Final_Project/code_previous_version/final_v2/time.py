import timeit


def f():
    #binary
    from simulation_data import N
    from constants import G, pi, msun, au, yr
    from physics import update
    from output import record
    import numpy as np


    # initial setup
    step     = 0             #start from 0th step
    dt       = 0.001*yr     #sec, time step
    time     = 0.0           #sec, start from t=0
    tmax     = 10*yr          #sec, max simulation time

    # output intervals to record trajectory
    noutput  = 200
    interval = (tmax/dt)/noutput

    # initialize the stars
    m1 = 1.0 * msun
    m2 = 2.0 * msun
    m3 = 1.0 * msun
    a           = 3*au
    a3          = 10*au
    period1     = (((4*pi**2)/(G*(m1+m2)))*a**3)**0.5
    period2     = (((4*pi**2)/(G*(m1+m2+m3)))*a3**3)**0.5
    force12     = (G*m1*m2)/(3*au)**2
    force13     = (G*m1*m3)/(8*au)**2
    force23     = (G*m3*m2)/(11*au)**2


    star1=np.array([1,        #id
                    m1,       #mass
                    -2.0*au,  #x
                    0,        #y
                    0,        #vx
                    (2.0*pi*m2/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2),
                    (-force13+force12)/m1, 
                    0         #ay
    ])

    star2=np.array([2,        #id
                    m2,       #mass
                    au,       #x
                    0,        #y
                    0,        #vx
                    (-2.0*pi*m1/(m1+m2)*a/period1)-(2.0*pi*(m3)/(m1+m2+m3)*a3/period2),
                    -(force12+force23)/m2, 
                    0         #ay
    ])

    star3=np.array([3,        #id
                    m3,       #mass
                    -10.0*au,  #x
                    0,        #y
                    0,        #vx
                    (2.0*pi*(m1+m2)/(m1+m2+m3)*a3/period2),
                    (force23+force13)/m3, 
                    0         #ay
    ])

    starnp=np.array([star1, star2, star3])
    first=True



    for n in np.arange(time, tmax, dt):

        starnp=update(starnp, dt)

        time=time+dt

        if (step%interval == 0):
            first=record(starnp, time, first)

        step = step+1


    print("Done!")

print("Using Python 5 times...")
t=timeit.timeit(f, number=5)
print("Using Python 5 times will take you[sec]:")
print(t)

