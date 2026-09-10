import numpy as np
import matplotlib.pyplot as plt

def fixed_point_iteration(g, p0, N):
    
    x = np.zeros((N, 1))
    
    x[0, 0] = p0
    
    for i in range(1, N):
        x[i, 0] = g(x[i-1, 0])
        
    return x