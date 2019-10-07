import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


f = open("some_output.out")
for i in range(1,6):
	f.readline()
i = f.readline()
x = []
y = []
for i in f:
	l = i.split(",")
	x.append(float(l[3]))
	y.append(float(l[4]))

xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

idx = z.argsort()
x1=[]
y1=[]
z1=[]
for i in range(len(x)):
	x1.append(x[idx[i]])
	y1.append(y[idx[i]])
	z1.append(z[idx[i]])

fig, ax = plt.subplots()
ax.scatter(x1, y1, c=z1, s=50, edgecolor='')
ax.annotate('Center', xy=(x1[-1],y1[-1]), xytext=(x1[-1],y1[-1]),
            arrowprops=dict(facecolor='black', shrink=0.05),
            )
plt.xlabel('L Label')
plt.ylabel('B Label')

plt.show()


