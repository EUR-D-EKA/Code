import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

# Fixing random state for reproducibility



f = open("../Output/hist_LMC_shift")
x = []
y = []
for i in f:
    if not i.startswith("1"):
        continue
    z = i.split()
    x.append(float(z[1]))
    y.append(float(z[3]))



fig, axs = plt.subplots()

# We can set the number of bins with the `bins` kwarg
axs.bar(x,y,align='center')

plt.ylim(0, 0.35)

plt.show()