import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

N, energy, estimate = np.loadtxt("efficiency.dat", unpack=True)
all_points = np.linspace(-1, 1, 1000)
all_exact = []

plt.figure(1)

energyplot   = plt.loglog(N, energy  , 'b-', label='energyplot')
estimateplot = plt.loglog(N, estimate, 'g-', label='estimateplot')
plt.grid(True)
plt.xlabel("norm")
plt.ylabel("N")
plt.legend()

plt.show()
#plt.figure(1).savefig("test.pdf")