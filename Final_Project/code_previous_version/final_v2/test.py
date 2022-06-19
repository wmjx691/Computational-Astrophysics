import numpy as np

x = np.array([[1,2,3],[4,5,6]])

a=7

np.savetxt('test.out', x) 
np.savetxt('test%d.out'%a, x[1],fmt='%1.4e') 
np.savetxt('test2.out', x, delimiter=',') 
np.savetxt('test3.out', x,newline='a') 
np.savetxt('test4.out', x,delimiter=',',newline='a') 
np.savetxt('test5.out', x,delimiter=',',header='abc') 
np.savetxt('test6.out', x,delimiter=',',footer='abc') 