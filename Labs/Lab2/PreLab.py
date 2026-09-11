import numpy as np
import matplotlib.pyplot as plt

def driver():
    # test functions
    f1 = lambda x: (10/(x+4))**(1/2)
    # fixed point is alpha1 = 1.4987....
    f2 = lambda x: 3+2*np.sin(x)
    #fixed point is alpha2 = 3.09...
    Nmax = 1000
    tol = 10**(-10)
    # test f1 '''
    p0 = 1.5
    [xstar,ier, count, pIterates] = fixed_point_iteration(f1, p0, tol, Nmax)
    print('the approximate fixed point is:',xstar)
    print('f1(xstar):',f1(xstar))
    print('count:', count)
    print('Error message reads:',ier)
    
def fixed_point_iteration(g, p0, tol, N):
    x = np.zeros((N,1))
    count  = 0
    while(count < N):
        count = count + 1
        p1 = g(p0)
        x[count-1, 0] = p1
        if(abs(p1 - p0) < tol):
            xstar = p1
            ier = 0
            return [xstar, ier, count, x]
        p0 = p1
    xstar = p1
    ier = 1
    return [xstar, ier, count, x]
        
driver()



