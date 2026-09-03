import numpy as np
import matplotlib.pyplot as plt

x = np.arange(1.92, 2.080, 0.001)

xPlot1 = [1, -18, 144, -672, 2016, -4032, 5376, -4608, 2304, -512]
xPlot2 = pow((x-2),9)

plt.plot(x, np.polyval(xPlot1, x), label='coefficents of (x-2)^9')
plt.plot(x, xPlot2, label='(x-2)^9')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

plt.show()

