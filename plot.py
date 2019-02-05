import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('/home/pratima/Landau_Ginzburg_SAM/output.txt',delimiter=' ')
t = data[:,0]
theta_z = data[:,1] * 180/np.pi

plt.figure(0)
plt.plot(t,theta_z)
plt.ylim([0,50])
plt.show()

print(theta_z[0])
print(theta_z[len(theta_z)-1])

