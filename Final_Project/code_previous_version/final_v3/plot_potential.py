#plot
from constants import *

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 

# initialize the stars
m1 = 1.0 * msun
m2 = 2.0 * msun
m3 = 1.0 * msun

r1 = -2.0
r2 = 1.0
r3 = -10.0

def V(x, y):
    return (G*m1/((x-r1)**2+y**2)**0.5+G*m2/((x-r2)**2+y**2)**0.5+G*m3/((x-r3)**2+y**2)**0.5)/au

n=1000
x1 = np.linspace(-15, 10, n)
y1 = np.linspace(-10, 10, n)

# x2=-x2[::-1]
# y2=-y2[::-1]

X, Y = np.meshgrid(x1, y1)

fig = plt.figure()

plt.title("The Gravitational Potential of 3-Body System")
 
plt.xlabel("X [AU]")
plt.ylabel("Y [AU]")

# contour=plt.contourf(X, Y, V(X, Y),levels=[6e12,6.8733e12,1e13,1.8205e13,2e13],norm = LogNorm())

cmap = mpl.colors.ListedColormap(['lightsteelblue',
                                'cornflowerblue',
                                'royalblue',
                                'navy'])
cmap.set_over('white')
cmap.set_under('white')

level=[3e12,6.8733e12,1e13,1.8205e13,3e13]
norm = mpl.colors.BoundaryNorm(level, cmap.N)
                            
plt.contourf(X, Y, V(X, Y),levels=level,norm = LogNorm(),cmap = cmap)

fig.colorbar(
    mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
    boundaries=[0]+level+[2e15],
    extend='both',
    ticks=level,
    spacing='proportional',
    orientation='vertical',
    label='Gravitational Potential [cm^2/sec^2]',
)

plt.savefig('Potential.png')
plt.show()