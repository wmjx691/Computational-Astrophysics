from constants import G, pi, msun, mass, au, yr
from solar import noutput, dt

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation


#1.input

df1 = pd.read_csv('planet_1.dat')
df2 = pd.read_csv('planet_2.dat')
df3 = pd.read_csv('planet_3.dat')
df4 = pd.read_csv('planet_4.dat')
df5 = pd.read_csv('planet_5.dat')
df6 = pd.read_csv('planet_6.dat')
df7 = pd.read_csv('planet_7.dat')
df8 = pd.read_csv('planet_8.dat')
df9 = pd.read_csv('planet_9.dat')
df10 = pd.read_csv('planet_10.dat')

# posx1=np.asarray(df1['posx [JPL]'])
# posx2=np.asarray(df2['posx [JPL]'])
posx3=np.asarray(df3['posx [JPL]'])/au
# posx4=np.asarray(df4['posx [JPL]'])
posx5=np.asarray(df5['posx [JPL]'])/au
posx6=np.asarray(df6['posx [JPL]'])/au
posx7=np.asarray(df7['posx [JPL]'])/au
posx8=np.asarray(df8['posx [JPL]'])/au
# posx9=np.asarray(df9['posx [JPL]'])
posx10=np.asarray(df10['posx [JPL]'])/au

# posy1=np.asarray(df1['posy [JPL]'])
# posy2=np.asarray(df2['posy [JPL]'])
posy3=np.asarray(df3['posy [JPL]'])/au
# posy4=np.asarray(df4['posy [JPL]'])
posy5=np.asarray(df5['posy [JPL]'])/au
posy6=np.asarray(df6['posy [JPL]'])/au
posy7=np.asarray(df7['posy [JPL]'])/au
posy8=np.asarray(df8['posy [JPL]'])/au
# posy9=np.asarray(df9['posy [JPL]'])
posy10=np.asarray(df10['posy [JPL]'])/au


Interval = (len(posx3)/dt)/noutput


#2.parameter setting

fig = plt.figure(figsize=(10, 10), dpi=50)
ax  = fig.gca()

ax.set_xlabel('x[AU]', fontsize=14)
ax.set_ylabel('y[AU]', fontsize=14)

ax.set_xlim(( -30, 10))
ax.set_ylim(( -30, 10))

#frame number
N = noutput


# 1Mercury
# 2Venus
# 3Earth
# 4Mars
# 5Jupiter
# 6Saturn
# 7Uranus
# 8Neptune
# 9Pluto
# 10Sun
 
#3.plot!

Sundot,     = ax.plot(posx10[:], posx10[:],
                color='red', marker='o', markersize=15,
                linestyle='')

Earthline,  = ax.plot([], [],
                color='green', linestyle='-', linewidth=1)
Earthdot,   = ax.plot([], [],
                color='green', marker='o', markersize=3,
                linestyle='')

Jupiterline,= ax.plot([], [],
                color='orange', linestyle='-', linewidth=1)
Jupiterdot, = ax.plot([], [],
                color='orange', marker='o', markersize=10,
                linestyle='')

Saturnline,= ax.plot([], [],
                color='brown', linestyle='-', linewidth=1)
Saturndot, = ax.plot([], [],
                color='brown', marker='o', markersize=8,
                linestyle='')

Uranusline,= ax.plot([], [],
                color='cyan', linestyle='-', linewidth=1)
Uranusdot, = ax.plot([], [],
                color='cyan', marker='o', markersize=5,
                linestyle='')

Neptuneline,= ax.plot([], [],
                color='blue', linestyle='-', linewidth=1)
Neptunedot, = ax.plot([], [],
                color='blue', marker='o', markersize=5,
                linestyle='')

step=0

def update(Interval):

    Earthline.set_data(posx3[0:Interval], posy3[0:Interval])
    Earthdot.set_data(posx3[Interval], posy3[Interval])
    Jupiterline.set_data(posx5[0:Interval], posy5[0:Interval])
    Jupiterdot.set_data(posx5[Interval], posy5[Interval])
    Saturnline.set_data(posx6[0:Interval], posy6[0:Interval])
    Saturndot.set_data(posx6[Interval], posy6[Interval])
    Uranusline.set_data(posx7[0:Interval], posy7[0:Interval])
    Uranusdot.set_data(posx7[Interval], posy7[Interval])
    Neptuneline.set_data(posx8[0:Interval], posy8[0:Interval])
    Neptunedot.set_data(posx8[Interval], posy8[Interval])

    print(Interval)

    return Earthline, Earthdot, Jupiterline, Jupiterdot, Saturnline, Saturndot, Uranusline, Uranusdot, Neptuneline, Neptunedot

def init():
    Earthline.set_data(posx3[0], posy3[0])
    Earthdot.set_data(posx3[0], posy3[0])
    Jupiterline.set_data(posx5[0], posy5[0])
    Jupiterdot.set_data(posx5[0], posy5[0])
    Saturnline.set_data(posx6[0], posy6[0])
    Saturndot.set_data(posx6[0], posy6[0])
    Uranusline.set_data(posx7[0], posy7[0])
    Uranusdot.set_data(posx7[0], posy7[0])
    Neptuneline.set_data(posx8[0], posy8[0])
    Neptunedot.set_data(posx8[0], posy8[0])

    return Earthline, Earthdot, Jupiterline, Jupiterdot, Saturnline, Saturndot, Uranusline, Uranusdot, Neptuneline, Neptunedot


ani = animation.FuncAnimation(fig=fig, func=update, frames=N, init_func=init, interval=365/N, blit=True, repeat=True)

plt.show()

#pillow is the writer of gif
ani.save('Solar System without stellite.gif', writer='pillow', fps=80)

print("DONE!!")