#simulation_data
import numpy as np
from constants import *

global N
N=3  #number of stars

class Simulation_data:
    def __init__(self, id, mass, x, y, vx, vy, ax, ay):
        self.id    =id
        self.mass  =mass
        self.x     =x
        self.y     =y
        self.vx    =vx
        self.vy    =vy
        self.ax    =ax
        self.ay    =ay

        return

class Stars(Simulation_data):
    pass

#     def __init__(self):
#         self=np.arange(1, N, dtype=Simulation_data)
    

# if __name__=='__main__':

