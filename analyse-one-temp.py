import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('log-300.txt', skip_header=1)
data = np.genfromtxt('output.300')

print 'total interaction energy  value: {}'.format(np.sum(data[:,2]))
print data.shape
tot_edges = np.where(data[:,2]>1)
print tot_edges

# print 'mean eta value: {}'.format(np.mean(data[1:,3]))
# print 'mean eta^2 value: {}'.format(np.mean(data[1:,4]))
