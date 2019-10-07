import matplotlib.pyplot as plt
from math import sqrt , log , cos ,sin,pi,log10,e
import scipy.integrate as integrate
from numpy import inf
from scipy.special import erf
import copy
from Code.Func_Class.Point_Class import Point

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
	k = 0.68
	r = point.r()
	Lambda = r/1.6/k
	a = 12
	v0 = 73
	
	#sigma2=integrate.quad(lambda x: 4*pi/3/v(point.r())*(x**4)*Distribution_function(0.5*x**2+FH), 0, 1000)
	#sigma = sqrt(sigma2[0])	
	sigma = isotropic_v_dispersion(a,v0,r)
	X = point.v()/sqrt(2*sigma)
	
	density = 5329/(2*pi)*((3+r**2/144)/144/(1+r**2/144)**2)

	F = 4*pi*point.m**2*log(Lambda)*density/point.v()**2*(erf(X)-2*X/pi**(1/2)*e**(-X**2))
	return (F*point.vx/point.v(), F*point.vy/point.v(), F*point.vz/point.v(), F)

def acceration(point):
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

mass = []
f = open('mass_file.txt')
f2 = open('center_fileWithoutDF','w')
f3 = open('dynamical_friction','w')
for i in f:
	mass.append(float(i))

l ,b ,r = 218 , 53.5 , 28.6
x ,y ,z = -15.704277369191800,      6.221783663213186,     15.055919236781001
vx ,vy , vz = 154.350088744764065,    166.787125980052366,   -106.603962915386660
point = Point(x,y,z,vx,vy,vz,mass[0])
X = [point.x]
Y = [point.y]
Z = [point.z]
dynamical_frictions = [0]

for i in range(1,7680):
	a = acceration(point)
	fake_point = copy.copy(point)
	fake_point.probe_advance(a[0],a[1],a[2], 0.00051398136723478,mass[i])
	a = acceration(fake_point)
	point.advance(a[0],a[1],a[2], 0.00051398136723478,mass[i])
	X.append(point.x)
	Y.append(point.y) 
	Z.append(point.z)
	dynamical_frictions.append(a[3])

for i in range(len(X)):
	f2.write("{} {} {}\n".format(X[i],Y[i],Z[i]))
	f3.write("{}\n".format(dynamical_frictions[i]))
'''
plt.ion()
animated_plot = plt.plot(X, Y, 'ro')[0]

for i in range(len(X)):
    animated_plot.set_xdata(X[0:i])
    animated_plot.set_ydata(Y[0:i])
    plt.draw()
    plt.pause(0.000000001)

for i in range(len(X)):
	print(X[i]," ",Y[i]," ",Z[i])
print("the final speed ",point.vx," ", point.vy," ",point.vz)
'''
