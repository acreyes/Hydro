import scipy as sp 
import numpy as np 
import matplotlib.pyplot as plt
import pdb
import os

data = np.loadtxt('data.dat')
X = data[:, 0]
rho = data[:,1]

plt.plot(X, rho)
plt.ylim(-.1, 1.1,)
plt.show()
