import numpy as np
import numpy.linalg as la

def eval_legendre(n, x):
    q = np.zeros(n+1)
    q[0] = 1
    
    if n >= 0:
        q[1] = x
        
    for i in range(1, n):
        q[i+1] = ((2*i + 1)*x*q[i] - i*q[i-1])/(i + 1)
        
    return q

