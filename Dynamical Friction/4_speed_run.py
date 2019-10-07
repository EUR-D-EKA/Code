import numpy as np
from scipy.stats import gaussian_kde
from math import sqrt
from timeit import default_timer as timer


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
		return (self.m*self.vx, self.m*self.vy, self.m*self.vz)
	def tp(self):
		p1 = self.p()
		return sqrt(p1[0]**2+p1[1]**2+p1[2]**2)
	def RKE(self,o):
		return 0.5 * self.m*((self.vx-o.vx)**2+(self.vy-o.vy)**2+(self.vz-o.vz)**2)
	def distance(self,o):
		return sqrt((self.x-o.x)**2+(self.y-o.y)**2+(self.z-o.z)**2)
	def bound(self,o):
		return (self.RKE(o)-self.m*o.m/self.distance(o)) < 0

def Center(points,x,y,z,r,m):
	p = [0,0,0]
	tm = 0
	fake = Point(x,y,z,0,0,0,0)
	for i in range(len(points)):
		if(Point.distance(points[i],fake)<=r):
			p[0] += points[i].p()[0]
			p[1] += points[i].p()[1]
			p[2] += points[i].p()[2]
			tm += points[i].m
		else:
			break
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
	save = 5
	for i in range(1,len(points)):
		tm += points[i].m
		if Point.bound(points[i],center):
			save = 5
		elif save == 0:
			x = tm
			for j in range(i,0,-1):
				x -= points[j].m
				if x <= tm*0.98:
					return Point.distance(points[j],center)
		else:
			save -=1
	x = tm
	for j in range(len(points)-1,-1,-1):
		x -= points[j].m
		if x <= tm*0.98:
			return Point.distance(points[j],center)

def data(file):
	f = open(file)
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
	point = []
	d =[]
	for i in range(len(x)):
		j = idx[i]
		d.append(density[j])
		point.append(Point(x[j],y[j],z[j],vx[j],vy[j],vz[j],m[j]))

	return point


m=600
r=13.4
totolm=[]
time = []
rad = []
center_ = []
for i in range(5,7764,6):
	start = timer()
	points = data(str(i))
	fake = Point(points[-1].x,points[-1].y,points[-1].z,0,0,0,0)
	points.sort(key = lambda p: Point.distance(p,fake))
	center = Center(points,points[0].x,points[0].y,points[0].z,r,m)
	r = radius(points,center)
	m = mass(points,center,r)
	center_.append((center.x,center.y,center.z))
	totolm.append(m)
	rad.append(r)
	time.append(i)
	end = timer()
	print("in progress!! ",i," time: ",end - start)

f = open("mass_fileSpeed","w")
f2 = open("center_fileSpeed","w")
f3 = open("radius_fileSpeed","w")

print('Center is ',center.x,' ',center.y,' ',center.z,' with radius ',r)
for i in range(len(time)):
	f.write(str(totolm[i])+"\n")
	f2.write("{} {} {}\n".format(center_[i][0],center_[i][1],center_[i][2]))
	f3.write(str(rad[i])+"\n")
	print("In timestep {} the center has mass is {}".format(time[i],totolm[i]))
