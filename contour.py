import matplotlib
import time
import math as math
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

traj = np.genfromtxt('/home/pratima/Landau_Ginzburg_SAM/traj1.txt')
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

size = 20
npts = len(traj)
x = np.linspace(0,19,20)
y = np.linspace(0,19,20)
X, Y = np.meshgrid(x, y)
z = np.zeros(size*size)
Z = np.zeros((size,size))

def assign_z(i):
    for k in range(0,400):
        z[k] = traj[i + k]
    return

plt.figure()
plt.ion()
plt.show()

for i in range(0,len(traj),400):
    assign_z(i)
    Z = np.reshape(z,(size,size))

    CS = plt.contour(X, Y, Z)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('y(r) at a specific time step')
#    plt.colorbar(CS, shrink=0.8, extend='both')
#    time.sleep(0.001)
    plt.pause(0.0001)


