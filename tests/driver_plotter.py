import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

points, approximate, exact = np.loadtxt("output.dat", unpack=True)
all_points = np.linspace(-1, 1, 1000)
all_exact = []

for i in range(0, 1000):
	x = all_points[i]
	all_exact.append(maths.sin(maths.pi * x))

plt.figure(1)

plt.plot(points,     approximate, 'b-')
plt.plot(all_points, all_exact,   'g-')
#plt.plot(points,     exact,       'r-')
plt.grid(True)
plt.xlabel("u")
plt.ylabel("x")

plt.show()
#plt.figure(1).savefig("test.pdf")