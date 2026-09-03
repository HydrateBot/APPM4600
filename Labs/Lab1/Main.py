import numpy as np
import matplotlib.pyplot as plt

x = [1,2,3]
y = np.array([1,2,3])

X = np.linspace(0, 2 * np.pi, 100)
Ya = np.sin(X)
Yb = np.cos(X)

plt.plot(X, Ya)
plt.plot(X, Yb)
plt.xlabel('x')
plt.ylabel('y')

plt.show()

