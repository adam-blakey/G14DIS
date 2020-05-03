import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math as maths

dof, error, indicator = np.loadtxt("../data/convergence.dat", unpack=True)

plt.figure(1)

plt.semilogy(dof,  error,                'bo-',                       label="Energy error")
plt.semilogy(dof,  indicator,            'go-',                       label="Error estimator")
plt.grid(True)
plt.xlabel("DoF")
plt.ylabel("error")
plt.title("Test Problem 1 with hp-Adaptivity Error")
plt.legend(loc="upper right")
print(indicator/error)
plt.show()