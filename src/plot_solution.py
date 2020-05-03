import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

points, approximate, exact = np.loadtxt("../data/solution.dat", unpack=True)

plt.figure(1)

plt.plot(points, approximate, 'b-', label="Approximation")
plt.plot(points, exact,       'g-', label="Exact")
plt.grid(True)
plt.xlabel("x")
plt.ylabel("u")
plt.title("Test Problem 1 with 1 p-Adaptive Step")
plt.legend(loc="lower left")

plt.show()