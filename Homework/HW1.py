import numpy as np
import matplotlib.pyplot as plt

x = np.arange(1.92, 2.080, 0.001)

xPlot = pow((x-9),9)

plt.plot(x, xPlot)
plt.xlabel('x')
plt.ylabel('y')


plt.show()

    