"""
 This script code compares the use of Lagrange and Newton Divided 
 Differences for interpolation. 
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
import matplotlib.pyplot as plt

def driver():


#    f = lambda x: np.exp(x)
    f = lambda x: 1/(1+(x**2))


    N = [5, 10, 15, 20] # number of interpolation nodes
    ''' interval'''
    a = -5
    b = 5
   
    for i, N in enumerate(N):
        ''' create equispaced interpolation nodes'''
        j = np.arange(N+1)
        xint = (a+b)/2 - (b-a)/2 * np.cos((2*j+1)*np.pi/(2*(N+1)))
        
        ''' create interpolation data'''
        yint = f(xint)
        
        ''' create points for evaluating the Lagrange interpolating polynomial'''
        Neval = 10000
        xeval = np.linspace(a,b,Neval+1)
        yeval_l= np.zeros(Neval+1) # lagrange
        
        '''Initialize and populate the first columns of the 
            divided difference matrix. We will pass the x vector'''
        y = np.zeros( (N+1, N+1) )
            
        
        ''' evaluate lagrange poly '''
        for kk in range(Neval+1):
            yeval_l[kk] = eval_lagrange(xeval[kk],xint,yint,N)
    
        ''' create vector with exact values'''
        fex = f(xeval)
            
    
        plt.figure()    
        plt.plot(xeval,fex,'ro-')
        plt.plot(xeval,yeval_l,'bs--') 
        plt.legend()
    
        plt.figure() 
        err_l = abs(yeval_l-fex)
        plt.semilogy(xeval,err_l,'ro--',label='lagrange')
        plt.legend()
        plt.show()

def eval_lagrange(xeval,xint,yint,N):

    lj = np.ones(N+1)
    
    for count in range(N+1):
       for jj in range(N+1):
           if (jj != count):
              lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    yeval = 0.
    
    for jj in range(N+1):
       yeval = yeval + yint[jj]*lj[jj]
  
    return(yeval)
  
    



       

driver()        
