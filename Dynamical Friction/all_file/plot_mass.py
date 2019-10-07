import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D


f = open("mass_file_10^9")
x=[]
t = 0
for line in f:
	x.append(float(line))


print(len(x),len(x))
plt.figure()
plt.plot(x)
plt.xlabel('Time Step')
plt.ylabel('Mass(structure unit)')
plt.title('Mass as function of time')
plt.show()