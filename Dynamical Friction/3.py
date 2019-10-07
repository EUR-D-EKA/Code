import matplotlib.pyplot as plt
from math import sqrt , log , cos ,sin,pi,log10


class Point(object):
	def __init__(self,x0,y0,z0,vx0,vy0,vz0,m0):
		self.x = x0
		self.y = y0
		self.z = z0
		self.vx = vx0
		self.vy = vy0
		self.vz = vz0
		self.m = m0  
	def advance(self,ax,ay,az,time,m0):
		self.x += self.vx*time + 0.5*ax*time**2
		self.y += self.vy*time + 0.5*ay*time**2
		self.z += self.vz*time + 0.5*az*time**2
		self.vx += ax*time
		self.vy += ay*time
		self.vz += az*time
		self.m = m0
	def r(self):
		return sqrt(self.x**2+self.y**2+self.z**2)
	def v(self):
		return sqrt(self.vx**2+self.vy**2+self.vz**2)
	def R(self):
		return sqrt(self.x**2+self.y**2)


def dynamical_friction(point):
	rho = 0
	k = 0.68
	Lambda = point.r/1.6/k
	sigma = 0
	X = point.v/sqrt(2*sigma)

	F = 4*pi*point.m**2*log(Lambda)*rho/point.v**2*(erf(X)-2*X/pi**(1/2)*exp(-X**2))
	return [F*point.vx/point.v, F*point.vy/point.v, F*point.vz/point.v]

def acceration(point):
	#return list that contain acceration in each axis
	FB = 449865.8882307301/((point.r()+0.7)**2)
	FDR = 152954.4019984482/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*point.R()
	FDZ = 152954.4019984482/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*(6.5+(point.z**2+0.0676)**0.5)/(point.z**2+0.0676)**0.5*point.z
	FH = 5329/(1+point.r()**2/144)*point.r()*2/144
	a = (FB + FH)
	ax = a * point.x/point.r() +FDR*point.x/point.R()
	ay = a * point.y/point.r() +FDR*point.y/point.R()
	az = a * point.z/point.r() +FDZ
	return [-ax,-ay,-az]

mass = []
f = open('mass_file.txt')
f2 = open('center_file2','w')
for i in f:
	mass.append(float(i))

l ,b ,r = 218 , 53.5 , 28.6
x ,y ,z = -15.325940787759894,5.125169260723238,15.681186049401791 
print(x,y,z)
vx ,vy , vz = -156, 79 ,107

point = Point(x,y,z,vx,vy,vz,mass[0])
X = [point.x]
Y = [point.y]
Z = [point.z]

for i in range(1,10):
	a = acceration(point)
	point.advance(a[0],a[1],a[2],0.00051398136723478,mass[i])
	X.append(point.x)
	Y.append(point.y)
	Z.append(point.z)

for i in range(len(X)):
	f2.write("{} {} {}\n".format(X[i],Y[i],Z[i]))

plt.ion()
animated_plot = plt.plot(X, Y, 'ro')[0]

for i in range(len(X)):
    animated_plot.set_xdata(X[0:i])
    animated_plot.set_ydata(Y[0:i])
    plt.draw()
    plt.pause(0.000000001)
for i in range(len(X)):
	print(X[i]," ",Y[i]," ",Z[i])