import matplotlib.pyplot as plt
import numpy as np

#I think this is all the prelab is

x0 = 1
x1 = 2


def line(x0, x1, fx, alpha):
    y0 = fx(x0)
    y1 = fx(x1)
    
    m = (y1 - y0) / (x1 - x0)
    
    y = y0 + m * (alpha - x0)
    
    return y