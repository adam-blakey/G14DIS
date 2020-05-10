import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

left, right, p = np.loadtxt("../data/mesh.dat", unpack=True)
points, approximate, exact = np.loadtxt("../data/solution.dat", unpack=True)

plt.figure(1)

plt.plot([left,  right], [p,     p],     'b')
plt.plot([left,  left],  [p-0.1, p+0.1], 'g')
plt.plot([right, right], [p-0.1, p+0.1], 'g')
#plt.plot(points,     approximate, 'r-')
#plt.plot(points, exact,   'y-')
plt.axis([0, 1, 0, max(p)+1])
plt.grid(True)
plt.xlabel("x")
plt.ylabel("p")
plt.title("Mesh Plot")


plt.show()