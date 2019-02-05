import pylab as plt
import numpy as np

X = np.linspace(0,19,20)
print X
Y = X**2 + np.random.random(X.shape)

plt.ion()
graph = plt.plot(X,Y)[0]
i = 0

while  i < 100:
    Y = X**2 + np.random.random(X.shape)
    graph.set_ydata(Y)
    plt.draw()
    plt.pause(0.01)
    i = i + 1
