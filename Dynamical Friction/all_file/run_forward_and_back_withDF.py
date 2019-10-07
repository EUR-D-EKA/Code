import matplotlib.pyplot as plt
from math import sqrt , log , cos ,sin,pi,log10,e
import scipy.integrate as integrate
from numpy import inf
from scipy.special import erf
import copy

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
	def probe_advance(self,ax,ay,az,time,m0):
		x = self.x + self.vx*time + 0.5*ax*time**2 
		y = self.y + self.vy*time + 0.5*ay*time**2 
		z = self.z + self.vz*time + 0.5*az*time**2
		self.x = (self.x + x)/2
		self.y = (self.y + y)/2
		self.z = (self.z + z)/2
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
'''
def Distribution_function(x):
	v = 73
	if x/v**2 > 700:
		return pi**(5/2)/sqrt(8)/144/v*inf*(2*sqrt(2)*inf*erf((sqrt(2)*sqrt(x))/v)+erf(sqrt(x)/v))
	f = pi**(5/2)/sqrt(8)/144/v*e**(x/v**2)*(2*sqrt(2)*e**(x/v**2)*erf((sqrt(2)*sqrt(x))/v)+erf(sqrt(x)/v))

	return f

def v(r):
	v = 73
	return v**2/4/pi/144*(3+r**2/144)/(1+r**2/144)**2
'''
def isotropic_v_dispersion(a,v0,r):
	return (v0**2)/2*(4*r*(a+r)**2*log(1+r/a)-a*r*(5*a+4*r))/(a**2*(2*a+r))

def anisotropic_v_dispersion(a,v0,r):
	return (v0**2)/2*(3*a+2*r)/(2*a+r)


def dynamical_friction(point):
	k = 1.428
	r = point.r()
	Lambda = r/1.6/k
	a = 12
	v0 = 73
	
	#sigma2=integrate.quad(lambda x: 4*pi/3/v(point.r())*(x**4)*Distribution_function(0.5*x**2+FH), 0, 1000)
	#sigma = sqrt(sigma2[0])	
	sigma = isotropic_v_dispersion(a,v0,r)
	X = point.v()/sqrt(2*sigma)
	
	density = 5329/(2*pi)*((3+r**2/144)/144/(1+r**2/144)**2)

	F = -4*pi*point.m**2*log(Lambda)*density/point.v()**2*(erf(X)-2*X/pi**(1/2)*e**(-X**2))
	return (F*point.vx/point.v(), F*point.vy/point.v(), F*point.vz/point.v(), abs(F))

def accelerationWithoutDF(point):
	#return list that contain acceration in each axis
	FB = 1.52954402e5/((point.r()+0.7)**2)
	FDR = 4.45865888e5/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*point.R()
	FDZ = 4.45865888e5/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*(6.5+(point.z**2+0.0676)**0.5)/(point.z**2+0.0676)**0.5*point.z
	FH = 5329/(1+point.r()**2/144)*point.r()*2/144
	a = (FB + FH)
	DF = dynamical_friction(point)
	ax = a * point.x/point.r() +FDR*point.x/point.R()
	ay = a * point.y/point.r() +FDR*point.y/point.R()
	az = a * point.z/point.r() +FDZ 
	return (-ax,-ay,-az, DF[3])

def accelerationWithDF(point):
	#return list that contain acceration in each axis
	FB = 1.52954402e5/((point.r()+0.7)**2)
	FDR = 4.45865888e5/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*point.R()
	FDZ = 4.45865888e5/((point.R()**2+(6.5+(point.z**2+0.0676)**0.5)**2)**(3/2))*(6.5+(point.z**2+0.0676)**0.5)/(point.z**2+0.0676)**0.5*point.z
	FH = 5329/(1+point.r()**2/144)*point.r()*2/144
	a = (FB + FH)
	DF = dynamical_friction(point)
	ax = a * point.x/point.r() +FDR*point.x/point.R()
	ay = a * point.y/point.r() +FDR*point.y/point.R()
	az = a * point.z/point.r() +FDZ 
	return (-ax+DF[0]/point.m,-ay+DF[1]/point.m,-az+DF[2]/point.m, DF[3])

mass = []
f = open('mass_file_10^10')
f2 = open('center_file_run_10^10','w')
f3 = open('dynamical_friction_10^10','w')

for i in f:
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))
	mass.append(float(i))

l ,b ,r = 218 , 53.5 , 28.6
x ,y ,z = -21.40558512, -10.47359097, 22.99030621
vx ,vy , vz = 156, -79, -107
point = Point(x,y,z,vx,vy,vz,mass[0])
X = [point.x]
Y = [point.y]
Z = [point.z]
dynamical_frictions = [0]

for i in range(1,77670):
	a = accelerationWithoutDF(point)
	fake_point = copy.copy(point)
	fake_point.probe_advance(a[0],a[1],a[2], 5.085429256646e-5,mass[i])
	a = accelerationWithoutDF(fake_point)
	point.advance(a[0],a[1],a[2], 5.085429256646e-5,mass[i])
	X.append(point.x)
	Y.append(point.y) 
	Z.append(point.z)
	dynamical_frictions.append(a[3])

point.m = mass[0]
point.vx = -point.vx
point.vy = -point.vy
point.vz = -point.vz
for i in range(1,77670):
	a = accelerationWithDF(point)
	fake_point = copy.copy(point)
	fake_point.probe_advance(a[0],a[1],a[2], 5.085429256646e-5,mass[i])
	a = accelerationWithDF(fake_point)
	point.advance(a[0],a[1],a[2], 5.085429256646e-5,mass[i])
	X.append(point.x)
	Y.append(point.y) 
	Z.append(point.z)
	dynamical_frictions.append(a[3])
print("center",point.x,point.y,point.z)
print("speed",point.vx,point.vy,point.vz)
for i in range(len(X)):
	f2.write("{} {} {}\n".format(X[i],Y[i],Z[i]))
	f3.write("{}\n".format(dynamical_frictions[i]))