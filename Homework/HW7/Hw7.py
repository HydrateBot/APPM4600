import numpy as np
from numpy.linalg import inv 
from numpy.linalg import norm 
import matplotlib.pyplot as plt

def f(x):
    return 1 / (1 + (10*x)**2)

def interpolationCoeff(x,y):
    N = len(x)
    V = np.vander(x, N)
    c = np.linalg.inv(V).dot(y)
    
    return c

def polynomialEval(c, x):
    p = np.zeros_like(x, dtype = float)
    for i in c:
        p = p*x + i
    return p

x_grid = np.linspace(-1, 1, 1001)
y_grid = f(x_grid)

for N in range(2,20):
    x = np.linspace(-1, 1, N)
    y = f(x)
    c = interpolationCoeff(x, y)
    p = polynomialEval(c, x_grid)
    print(f"N = {N}")
    print(f"Maximum p(x) = {np.max(p):.6f}")
    print(f"Minimum p(x) = {np.min(p):.6f}")
    print()

    if N in [2, 5, 10, 15, 16, 17, 18, 19]:
        plt.plot(x_grid, y_grid, label='f(x)')
        plt.plot(x_grid, p, label=f'P_{N-1}(x)')
        plt.plot(x, y, 'o', label='Interpolation Points')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.grid(True)
        plt.show()