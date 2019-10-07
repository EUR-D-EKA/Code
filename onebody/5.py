import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D


f = open("center_fileWithoutDF")
f2 = open("center_file")
x=[]
y=[]
z=[]
x1=[]
y1=[]
z1=[]
t = 0
for line in f:
	position = line.split()
	position1 = f2.readline().split()
	t+=1
	x.append(float(position[0]))
	y.append(float(position[1]))
	z.append(float(position[2]))
	x1.append(float(position1[0]))
	y1.append(float(position1[1]))
	z1.append(float(position1[2]))

print(len(x1),len(x))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.title('Red line is MW@home, blue line is mine')
ax.scatter(x,y,z,c='b', edgecolor='',s=1)
ax.scatter(x1,y1,z1,c='r', edgecolor='',s=1)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()
