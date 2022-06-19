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
    def __init__(self, idindex, mass, x, y, z, vx, vy, vz, rsoi):
        self.idindex = idindex
        self.mass    = mass
        self.x       = x
        self.y       = y
        self.z       = z
        self.vx      = vx
        self.vy      = vy
        self.vz      = vz
        self.rsoi    = rsoi
        
        return

class Stars(starinfo):
    pass

class speedupinfomation:
    def __init__(self, soi_ox, idindex, speed_ox, direction, dt_ratio, shortest=0):
        self.soi_ox     = soi_ox
        self.idindex    = idindex
        self.speed_ox   = speed_ox
        self.direction  = direction
        self.dt_ratio   = dt_ratio
        self.shortest   = shortest

        return


#     def __init__(self):
#         self=np.arange(1, N, dtype=Simulation_data)
    

# if __name__=='__main__':

    

