#record
from simulation_data import N

import numpy as np
from astropy.time import Time
import pandas as pd

def record(data, time, first, planet):

    if planet:
        time_str = np.empty(len(time), dtype='object')
        
        for i in range(len(time)):
            temp        = Time('%s' %time[i], format='jd')
            time_str[i] = temp.to_value('iso', subfmt='date') 
            
    
        for i in np.arange(0, N):
            nptag = np.full((len(time)), data[i].idindex)
            npmass = np.full((len(time)), data[i].mass)
            
        
            df = pd.DataFrame({'Time [JD]'  : time,
                               'Time [UTC]' : time_str,
                               'Tag'        : nptag,
                               'Mass [kg]'  : npmass,
                               'posx [JPL]' : data[i].x,
                               'posy [JPL]' : data[i].y,
                               'posz [JPL]' : data[i].z,
                               'velx [JPL]' : data[i].vx,
                               'vely [JPL]' : data[i].vy,
                               'velz [JPL]' : data[i].vz,
                               })   

            #drop lanuch day's data
            df=df.drop([0])                     
        
            with open('planet_%d.dat' %(i+1), 'w') as f:
                df.to_csv(f, index=False)
                f.close()
        
        print("record each planet in solar system...")

        return
    
    else:
        temp     = Time('%s' %time, format='jd')
        time_str = temp.to_value('iso', subfmt='date') 

        if first:
            with open('stell.dat', 'w') as f:
                f.write("Time [JD]")
                f.write(',') 
                f.write("Time [UTC]")
                f.write(',') 
                f.write("Tag")
                f.write(',') 
                f.write("Mass [kg]")
                f.write(',')
                f.write("posx [code]")
                f.write(',') 
                f.write("posy [code]")
                f.write(',') 
                f.write("posz [code]")
                f.write(',') 
                f.write("velx [code]")
                f.write(',') 
                f.write("vely [code]")
                f.write(',') 
                f.write("velz [code]")
                f.write(',') 
                f.write("accx [code]")
                f.write(',') 
                f.write("accy [code]")
                f.write(',') 
                f.write("accz [code]")
                f.write('\n')
                f.close()

            first = False
            return first

        else:
            with open('stell.dat', 'at') as f:
                f.write(str(time))
                f.write(',')
                f.write(str(time_str))
                f.write(',')
                f.write(str(data.idindex))
                f.write(',')
                f.write(str(data.mass))
                f.write(',')
                f.write(str(data.x))
                f.write(',')
                f.write(str(data.y))
                f.write(',')
                f.write(str(data.z))
                f.write(',')                
                f.write(str(data.vx))
                f.write(',')
                f.write(str(data.vy))
                f.write(',')
                f.write(str(data.vz))
                f.write(',')
                f.write(str(data.ax))
                f.write(',')
                f.write(str(data.ay))
                f.write(',')
                f.write(str(data.az))
                #np.savetxt(f, data, newline=' ', fmt='%1.8e')
                f.write('\n')
                f.close()

            return first