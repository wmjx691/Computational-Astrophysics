#simulation_data
import numpy as np
from constants import G, pi, msun, mass, au, yr

global N
N=10        #number of stars

class stell_info:
    #def __init__(self, idindex, mass, x, y, z, vx, vy, vz, ax, ay, az):
    def __init__(self, idindex, mass, x, y, z, vx, vy, vz, ax, ay, az):
        self.idindex = idindex
        self.mass    = mass
        self.x       = x
        self.y       = y
        self.z       = z
        self.vx      = vx
        self.vy      = vy
        self.vz      = vz
        self.ax      = ax
        self.ay      = ay
        self.az      = az

        return
        
class starinfo:
    #def __init__(self, idindex, mass, x, y, z, vx, vy, vz, ax, ay, az):
    def __init__(self, idindex, mass, x, y, z, vx, vy, vz):
        self.idindex = idindex
        self.mass    = mass
        self.x       = x
        self.y       = y
        self.z       = z
        self.vx      = vx
        self.vy      = vy
        self.vz      = vz
        
        return

class Stars(starinfo):
    pass



#     def __init__(self):
#         self=np.arange(1, N, dtype=Simulation_data)
    

if __name__=='__main__':
    earth = star_info(1)
    

