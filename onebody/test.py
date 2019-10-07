import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from math import sqrt

f = open("center_file")
f2 = open("center_fileWithoutDF")
f3 = open("dynamical_friction")
distance = []
time_step = []
DF = []
i = 0
for line in f2:
	position = line.split()
	real_position = f.readline().split()
	DF.append(float(f3.readline()))
	distance.append(sqrt((float(position[0])-float(real_position[0]))**2+(float(position[1])-float(real_position[1]))**2+(float(position[2])-float(real_position[2]))**2))
	time_step.append(i)
	i += 1
print(len(distance),len(time_step))

plt.hist(distance)
plt.xlabel('Difference')
plt.ylabel('Error Frequency')
plt.title('distance error')
plt.grid(True)

plt.figure()
plt.plot(time_step,distance)
plt.title('Difference (3.95Gyr with mass of 96)')
plt.xlabel('time_step')
plt.ylabel('Difference(kpc)')

plt.figure()
plt.plot(time_step,DF)
plt.title('DF for one body(3.95Gyr with mass of 96)')
plt.xlabel('time_step')
plt.ylabel('unit force')


plt.show()
