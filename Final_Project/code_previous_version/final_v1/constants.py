#constant
import math

global G, pi, msun, mass, au, yr

s    = 86400.                  #day->s
m    = 1000.                   #km->m
cm   = 100000.0                #km->cm
kg   = 1000.0                  #kg->g

G    = 6.67428e-11*(s**2)/(m**3) #(km^3*kg^-1*day^-2)
pi   = 4.0*math.atan(1.0)
msun = 1.989e30                #(kg)
mass = 10.e24                  #(kg)
au   = 1.49598e8               #(km)
yr   = 365.                    #(d)


    # def show_constants(self):
    #     if self == "G":
    #         print(self.G)
    #     elif self == "pi":
    #         print(self.pi)
    #     elif self == "msun":
    #         print(self.msun)
    #     elif self == "au":
    #         print(self.au)
    #     elif self == "yr":
    #         print(self.yr)
        
    #     return