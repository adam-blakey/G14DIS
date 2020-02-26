import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

x, l0, l1, l2, l3, l4, l5, l6, l7, l8, l9 = np.loadtxt("plotTest.dat", unpack=True)

plt.figure(1)

plt.plot(x, l0, 'k-')
plt.plot(x, l1, 'k-')
plt.plot(x, l2, 'k-')
plt.plot(x, l3, 'r-')
plt.plot(x, l4, 'g-')
plt.plot(x, l5, 'b-')
plt.grid(True)
plt.xlabel("l")
plt.ylabel("x")

plt.show()
#plt.figure(1).savefig("test.pdf")