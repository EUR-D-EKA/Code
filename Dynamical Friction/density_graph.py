import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde

f = open("0")
for i in range(1,6):
	f.readline()
i = f.readline()

x,y,z,vx,vy,vz,m = [],[],[],[],[],[],[]
for i in f:
	l = i.split(",")
	x.append(float(l[2]))
	y.append(float(l[3]))
	z.append(float(l[4]))
	vx.append(float(l[8]))
	vy.append(float(l[9]))
	vz.append(float(l[10]))
	m.append(float(l[11]))

xyz = np.vstack([x,y,z])
kde = gaussian_kde(xyz)
density = kde(xyz)

idx = density.argsort()
x1,y1,z1,vx1,vy1,vz1,m1,d = [],[],[],[],[],[],[],[]
for i in range(len(x)):
	j = idx[i]
	x1.append(x[j])
	y1.append(y[j])
	z1.append(z[j])
	vx1.append(vx[j])
	vy1.append(vy[j])
	vz1.append(vz[idx[i]])
	m1.append(m[j])
	d.append(density[j])

print("Center is",x1[-1],y1[-1],z1[-1])
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1, y1, z1, c=d,edgecolor='',s = 1)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()


