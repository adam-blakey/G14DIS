import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

points, approximate, exact = np.loadtxt("output.dat", unpack=True)

plt.figure(1)

plt.plot(points, approximate, 'g-')
plt.plot(points, exact,       'b-')
plt.grid(True)
plt.xlabel("u")
plt.ylabel("x")

plt.show()
#plt.figure(1).savefig("test.pdf")