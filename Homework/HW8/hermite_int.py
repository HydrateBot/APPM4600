"""
 This script compares the use of Lagrange and Hermite interpolation.
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



import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import math

def driver():

    f = lambda x: 1./(1.+x**2)
    fp = lambda x: -2*x/(1.+x**2)**2

    N_values = [5, 10, 15, 20]

    ''' interval '''
    a = -5
    b = 5

    for N in N_values:

        ''' create equispaced interpolation nodes '''
        j = np.arange(N+1)
        xint = (a+b)/2 - (b-a)/2 * np.cos((2*j+1)*np.pi/(2*(N+1)))

        ''' create interpolation data '''
        yint = np.zeros(N+1)
        ypint = np.zeros(N+1)

        for jj in range(N+1):
            yint[jj] = f(xint[jj])
            ypint[jj] = fp(xint[jj])

        ''' create evaluation points '''
        Neval = 1000
        xeval = np.linspace(a,b,Neval+1)
        yevalH = np.zeros(Neval+1)

        ''' evaluate Hermite interpolating polynomial '''
        for kk in range(Neval+1):
            yevalH[kk] = eval_hermite(
                xeval[kk],xint,yint,ypint,N
            )

        ''' create vector with exact values '''
        fex = np.zeros(Neval+1)

        for kk in range(Neval+1):
            fex[kk] = f(xeval[kk])

        ''' graph Hermite interpolation '''
        plt.figure()
        plt.plot(xeval,fex,'r-',label='Exact Function')
        plt.plot(xeval,yevalH,'c--',label='Hermite')
        plt.plot(xint,yint,'ko',label='Interpolation Nodes')

        plt.title(f'Hermite Interpolation N = {N}')
        plt.xlabel('x')
        plt.ylabel('f(x)')
        plt.legend()
        plt.grid()

        ''' calculate Hermite error '''
        errH = abs(yevalH-fex)

        ''' graph Hermite error '''
        plt.figure()
        plt.semilogy(xeval,errH,'c-',label='Hermite Error')

        plt.title(f'Hermite Interpolation Error N = {N}')
        plt.xlabel('x')
        plt.ylabel('Absolute Error')
        plt.legend()
        plt.grid()

    plt.show()


def eval_hermite(xeval,xint,yint,ypint,N):

    ''' Evaluate all Lagrange polynomials '''

    lj = np.ones(N+1)

    for count in range(N+1):
        for jj in range(N+1):
            if (jj != count):
                lj[count] = lj[count]*(xeval - xint[jj])/(xint[count]-xint[jj])

    ''' Construct the l_j'(x_j) '''

    lpj = np.zeros(N+1)

    for count in range(N+1):
        for jj in range(N+1):
            if (jj != count):
                lpj[count] = lpj[count] + 1./(xint[count] - xint[jj])

    yeval = 0.

    for jj in range(N+1):
        Qj = (1.-2.*(xeval-xint[jj])*lpj[jj])*lj[jj]**2
        Rj = (xeval-xint[jj])*lj[jj]**2

        yeval = yeval + yint[jj]*Qj + ypint[jj]*Rj

    return(yeval)


if __name__ == '__main__':
    driver()
        
