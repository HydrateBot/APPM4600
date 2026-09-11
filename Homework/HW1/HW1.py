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

# iii.) When using the coefficients, python is forced to complete
# each operation in the polynomial, which can lead to compounding floating
# point errors. When using the (x-2)^9 form, python is able to compute
# the value of (x-2)^9 using only two operations, (x-2) then raising it to the 9th power.
# which significantly reduces the floating point error. Thus the correct graph
# is the graph of (x-2)^9 as opposed to the graph of the coefficients.

x2 = [np.pi, 10**6]
delta = np.logspace(-16, 0, 17)
for c in x2:
    oExpr = np.cos(c + delta) - np.cos(c)
    sExpr = -2 * np.sin(c + delta / 2) * np.sin(delta / 2)
    myExpr = -delta * np.sin(c)
    diff = np.abs(oExpr - sExpr)
    diff2 = np.abs(myExpr - sExpr)
    #plt.plot(delta, diff)
    plt.plot(delta, diff2)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('delta (log scale)')
plt.ylabel('Absolute Difference')
plt.legend()
plt.show()





