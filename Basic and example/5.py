import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import gaussian_kde
from math import sqrt


plt.figure(1)
plt.plot([1,2,3],[2,8,4])
plt.axis([0, 10, 0, 10])

plt.savefig('foo1.png')

plt.figure(2)
plt.plot([1,2,3],[9,1,2])
plt.axis([0, 10, 0, 10])

plt.savefig('foo2.png')
plt.show()