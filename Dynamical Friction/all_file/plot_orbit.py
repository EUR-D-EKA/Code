import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D


f = open("center_file_run_10^10")
x=[]
y=[]
z=[]
for line in f:
	position = line.split()
	x.append(float(position[0]))
	y.append(float(position[1]))
	z.append(float(position[2]))

print(len(x))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,z,c='b', edgecolor='',s=1)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()