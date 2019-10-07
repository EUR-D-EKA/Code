import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from math import sqrt
from scipy.constants import G




f = open("some_output.out")
for i in range(1,6):
	f.readline()
i = f.readline()
x = []
y = []
z = []
m = []
for i in f:
	l = i.split(",")
	x.append(float(l[1]))
	y.append(float(l[2]))
	z.append(float(l[3]))
	m.append(float(l[10]))


xyz = np.vstack([x,y,z])
kde = gaussian_kde(xyz)
density = kde(xyz)

idx = density.argsort()
x1=[]
y1=[]
z1=[]

m1= []
d=[]
for i in range(len(x)):
	x1.append(x[idx[i]])
	y1.append(y[idx[i]])
	z1.append(z[idx[i]])
	m1.append(m[idx[i]])
	d.append(density[idx[i]])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1, y1, z1,c=d, edgecolor='')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
print('Center is ',x1[-1],' ',y1[-1],' ',z1[-1])
plt.show()

