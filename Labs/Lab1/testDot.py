"""
 This program is a warm up for coding. You get used to the coding 
format and practice some coding skills. 
"""

############################################# 
"""
Copyright (C) 2025  Adrianna M. Gillman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
############################################# 


import numpy as np
import numpy.linalg as la
import math

def driver():

     n = 2
     x = np.linspace(0,np.pi,n)
     matrix1 = [[2,2],
                [2,2]]
     matrix2 = [[3,3],
                [3,3]]

# this is a function handle.  You can use it to define 
# functions instead of using a subroutine like you 
# have to in a true low level language.     
     f = lambda x: 2*x
     g = lambda x: 0*x

     y = f(x)
     w = g(x)
     
     print(y)
     print(w)

# evaluate the dot product of y and w     
     dp = dotProduct(y,w,n)
     dp2 = dotProduct(matrix1, matrix2, 2)
     

# print the output
     print('the dot product is : ', dp)
     print(dp2)

     return
     
def dotProduct(x,y,n):
#   Computes the dot product of the n x 1 vectors x and y
     dp = 0.
     for j in range(n):
        dp = dp + x[j]*y[j]

     return dp  
     
driver()     

            
