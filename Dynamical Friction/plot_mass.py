import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D


f = open("center_fileSpeed")
x=[]
t = 0
for line in f:
	real_position = f.readline().split()
	distance.append(sqrt((float(position[0])**2+float(position[1])**2+float(position[2])**2))
	t+=1

print(len(x),len(x))
plt.figure()
plt.plot(t,x)
plt.xlabel('Gyr')
plt.ylabel('Mass')
plt.show()