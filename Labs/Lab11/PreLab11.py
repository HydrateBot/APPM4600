import numpy as np

def cTrap(a,b,f,N):
    h = (b - a) / (N-1)
    x = np.linspace(a, b, N)
    
    I = f(x[0]) + f(x[-1])
    
    for i in range(1, N-1):
        I += 2 * f(x[i])
    
    return (h/2) * I

def simpsons(a,b,f,N):
    h = (b - a) / (N-1)
    x = np.linspace(a, b, N)
    
    I = f(x[0]) + f(x[-1])
    
    for i in range(1, N-1):
        if i % 2 == 0:
            I += 2 * f(x[i])
        else:
            I += 4 * f(x[i])
    
    return (h/3) * I



