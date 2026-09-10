import numpy as np
import matplotlib.pyplot as plt

def fixed_point_iteration(g, p0, tol, N):
    
    count  = 0
    while(count < N):
        count = count + 1
        p1 = g(p0)
        if(abs(p1 - p0) < tol):
            xstar = p1
            ier = 0
            return [xstar, ier]
        p0 = p1
    xstar = p1
    ier = 1
    return [xstar, ier]
        
    return x