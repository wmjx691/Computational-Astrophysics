#record
from simulation_data import N
import numpy as np

def record(starnp, time, first):
    for i in np.arange(0, N):
        if first:
            with open('binary_%d.dat' %(i+1), 'w') as f:
                f.write("Time [sec]")
                f.write('         ,') 
                f.write("# Tag")
                f.write('         ,') 
                f.write("Mass [g]")
                f.write('     ,') 
                f.write("posx [code]")
                f.write('    ,') 
                f.write("posy [code]")
                f.write('   ,') 
                f.write("velx [code]")
                f.write('    ,') 
                f.write("vely [code]")
                f.write('    ,') 
                f.write("accx [code]")
                f.write('    ,') 
                f.write("accy [code]")
                f.write('\n')
                f.close()
        else:
            with open('binary_%d.dat' %(i+1), 'at') as f:
                f.write(str(time))
                f.write(' ')                
                np.savetxt(f, starnp[i], newline=' ', fmt='%1.8e')
                f.write('\n')
                f.close()
        
    first=False
    return first