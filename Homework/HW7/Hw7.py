import numpy as np
from numpy.linalg import inv 
from numpy.linalg import norm 
import matplotlib.pyplot as plt
from math import comb

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

def barycentricWeights(N):
    n = N - 1
    
    w = np.array([(-1)**j * comb(n, j) for j in range(N)], dtype = float)
    w = w / np.max(np.abs(w))
    return w

def barycentricEval(x, y, w, x_eval):
    N = len(x)
    p_eval = np.zeros_like(x_eval, dtype = float)
    for i, x_i in enumerate(x_eval):
        diff = x_i - x
        match = np.where(np.abs(diff) < 1e-14)[0]

        if(len(match) > 0):
            p_eval[i] = y[match[0]]
        else:
            terms = w / diff
            num = np.sum(terms * y)
            denom = np.sum(terms)
            p_eval[i] = num / denom
    return p_eval

def interpolateBC(N, x_val):
    x = np.linspace(-0.2, 0.2, N)
    y = f(x)
    w = barycentricWeights(N)
    p_eval = barycentricEval(x, y, w, x_val)
    return x, y, p_eval

def chebyshevWeights(x):
    N = len(x)
    w = np.zeros(N)
    
    for j in range(N):
        product = 1
        for i in range(N):
            if i != j:
                product *= (x[j] - x[i])
        w[j] = 1 / product
    return w

def interpolateChebyshev(N, x_val):
    j = np.arange(1, N+1)
    x = np.cos((2*j - 1) * np.pi / (2*N))
    y = f(x)
    w = chebyshevWeights(x)
    p_eval = barycentricEval(x, y, w, x_val)
    return x, y, p_eval

x_grid = np.linspace(-1, 1, 1001)
y_grid = f(x_grid)

#For problem 3

for N in [5, 15, 20, 50, 2000]:
    x, y, p_eval = interpolateChebyshev(N, x_grid)
    
    plt.plot(x_grid, y_grid, label='f(x)')
    plt.plot(x_grid, p_eval, label=f'P_{N-1}(x)')
    plt.plot(x, y, 'o', label='Interpolation Points')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.title(f'Chebyshev Interpolation with N={N}')
    plt.show()


#For problem 2

# for N in [5, 15, 20, 50]:
#     x, y, p_eval = interpolateBC(N, x_grid)
    
#     plt.plot(x_grid, y_grid, label='f(x)')
#     plt.plot(x_grid, p_eval, label=f'P_{N-1}(x)')
#     plt.plot(x, y, 'o', label='Interpolation Points')
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.legend()
#     plt.grid(True)
#     plt.title(f'Barycentric Interpolation with N={N}')
#     plt.show()


# For problem 1

# for N in range(2,20):
#     x = np.linspace(-1, 1, N)
#     y = f(x)
#     c = interpolationCoeff(x, y)
#     p = polynomialEval(c, x_grid)
#     print(f"N = {N}")
#     print(f"Maximum p(x) = {np.max(p):.6f}")
#     print(f"Minimum p(x) = {np.min(p):.6f}")
#     print()

#     if N in [2, 5, 10, 15, 16, 17, 18, 19]:
#         plt.plot(x_grid, y_grid, label='f(x)')
#         plt.plot(x_grid, p, label=f'P_{N-1}(x)')
#         plt.plot(x, y, 'o', label='Interpolation Points')
#         plt.xlabel('x')
#         plt.ylabel('y')
#         plt.legend()
#         plt.grid(True)
#         plt.show()