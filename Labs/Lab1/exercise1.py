import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 6, 3)
y = np.arange(3)

print(x[0])

print('the first three entries of x are:', x[0], x[1], x[2])

w = 10**(-np.linspace(1,10,10))

x2 = np.linspace(1, 10, w.size,dtype=int)
s = 3 * w
plt.semilogy(x2, w, label ='x2 versus w')
plt.semilogy(x2, s, label ='x2 versus s')
plt.xlabel('x2')
plt.ylabel('w')
plt.legend()
plt.show()


