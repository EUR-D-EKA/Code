import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from math import sqrt


class Point(object):
	def __init__(self,x0,y0,z0,vx0,vy0,vz0,m0):
		self.x = x0
		self.y = y0
		self.z = z0
		self.vx = vx0
		self.vy = vy0
		self.vz = vz0
		self.m = m0
	def p(self):
		return [self.m*self.vx, self.m*self.vy, self.m*self.vz]
	def tp(self):
		p1 = self.p()
		return sqrt(p1[0]**2+p1[1]**2+p1[2]**2)
	def RKE(self,o):
		return 0.5 * self.m*((self.vx-o.vx)**2+(self.vy-o.vy)**2+(self.vz-o.vz)**2)
	def distance(self,o):
		return sqrt((self.x-o.x)**2+(self.y-o.y)**2+(self.z-o.z)**2)
	def bound(self,o):
		return (self.RKE(o)-self.m*o.m/self.distance(o)) < 0

def Center(points,x,y,z,m):
	p = [0,0,0]
	tm = 0
	for i in range(len(points)):
		p[0] += points[i].p()[0]
		p[1] += points[i].p()[1]
		p[2] += points[i].p()[2]
		tm += points[i].m
	return Point(x,y,z,p[0]/tm,p[1]/tm,p[2]/tm,m)

def mass(point,center,r):
	j = 0
	tm = 0
	points = point
	for i in range(len(points)):
		if Point.distance(points[i],center) <= r :
			j = i
		else:
			break
	for i in range(j+1):
		tm += points[i].m
	return tm

def radius(point,center):
	tm = 0
	x =0
	points = point
	save = 20
	for i in range(1,len(points)):
		tm += points[i].m
		if Point.bound(points[i],center) :
			save = 20
		elif save == 0:
			break
		else:
			save -=1
	for i in range(len(points)):
		x +=  points[i].m
		if x >= tm*0.98:
			return Point.distance(points[i],center)

def data(file):
	f = open(file)
	for i in range(1,6):
		f.readline()
	i = f.readline()

	x,y,z,vx,vy,vz,m = [],[],[],[],[],[],[]
	for i in f:
		l = i.split(",")
		x.append(float(l[1]))
		y.append(float(l[2]))
		z.append(float(l[3]))
		vx.append(float(l[7]))
		vy.append(float(l[8]))
		vz.append(float(l[9]))
		m.append(float(l[10]))

	xyz = np.vstack([x,y,z])
	kde = gaussian_kde(xyz)
	density = kde(xyz)

	idx = density.argsort()
	x1,y1,z1,vx1,vy1,vz1,m1,d = [],[],[],[],[],[],[],[]
	for i in range(len(x)):
		x1.append(x[idx[i]])
		y1.append(y[idx[i]])
		z1.append(z[idx[i]])
		vx1.append(vx[idx[i]])
		vy1.append(vy[idx[i]])
		vz1.append(vz[idx[i]])
		m1.append(m[idx[i]])
		d.append(density[idx[i]])

	point = []
	for i in range(len(x)):
		point.append(Point(x1[i],y1[i],z1[i],vx1[i],vy1[i],vz1[i],m1[i]))
	re = [point,d,x1,y1,z1]
	return re



m=60
totolm=[60]
time = [0]
for i in range(7000,7001):
	data1 = data(str(i))
	center = Center(data1[0],data1[0][-1].x,data1[0][-1].y,data1[0][-1].z,m)
	points = data1[0]
	points.sort(key = lambda p: Point.distance(p,center))
	r = radius(points,center)
	m = mass(points,center,r)
	totolm.append(m)
	time.append(i+1)
	if i == 7000:
		point = data1[0]
		d = data1[1]
		x = data1[2]
		y = data1[3]
		z = data1[4]


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z,c=d, edgecolor='')
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
print('Center is ',center.x,' ',center.y,' ',center.y,' with radius ',r)
for i in range(len(time)):
	print("In timestep {} the center has mass is {}".format(time[i],totolm[i]))
plt.show()

