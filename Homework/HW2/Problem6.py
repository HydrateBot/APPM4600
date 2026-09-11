import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x - 4*np.sin(2*x) - 3

x = np.linspace(-2, 8, 1000)
y = f(x)


plt.plot(x, y, label='f(x) = x - 4*sin(2*x) - 3')
plt.xlabel('x')
plt.ylabel('y')
plt.grid()
plt.legend()

plt.show()