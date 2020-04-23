import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

dof, error, indicator = np.loadtxt("convergence.dat", unpack=True)

plt.figure(1)

plt.semilogy(dof, error,     'bo-')
plt.semilogy(dof, indicator, 'go-')
plt.grid(True)
plt.xlabel("DoF")
plt.ylabel("error")
print(indicator/error)
plt.show()