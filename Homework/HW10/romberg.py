"""
 This script tests Romberg integration for approximating integrals.
 The first column of the Romberg table is built with the composite
 trapezoidal rule and the remaining columns come from Richardson
 extrapolation.  The error in the diagonal entries is plotted as the
 number of panels is increased.
"""

#############################################
"""
Copyright (C) 2026  Adrianna M. Gillman

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
from scipy.integrate import quad


def driver():
    f = lambda x: 1/(1+x**2) 
    a = -5
    b = 5
    
    I_ex = 2*np.arctan(5)
    tol = 1e-4
    
    M2 = 2
    M4 = 24
    
    ntrap = int(np.ceil(np.sqrt((b-a)**3*M2/(12*tol))))
    nsimp = int(np.ceil(((b-a)**5*M4/(180*tol))**(1/4)))
    
    if nsimp % 2 != 0:
        nsimp += 1
    
    iTrap = CompTrap(a,b,ntrap,f)
    iSimp = CompSimp(a,b,nsimp,f)
    
    print('trapezoidal rule with n =', ntrap, 'panels gives I =', iTrap, 'with error =', abs(I_ex-iTrap), 'function evaluations =', ntrap+1)
    print('Simpson''s rule with n =', nsimp, 'panels gives I =', iSimp, 'with error =', abs(I_ex-iSimp), 'function evaluations =', nsimp+1)
    
    res1 = quad(f,a,b,epsabs=1e-6, epsrel=1e-6,full_output=1)
    res2 = quad(f,a,b,epsabs=1e-4, epsrel=1e-4,full_output=1)
    
    print('quad with tolerance 1e-6 gives I =', res1[0], 'with function evaluations =', res1[2]['neval'])
    print('quad with tolerance 1e-4 gives I =', res2[0], 'with function evaluations =', res2[2]['neval'])

    # f = lambda x: np.exp(x) * np.cos(x)
    # a = 0
    # b = np.pi

    # # exact integral
    # I_ex = -12.070346316389632
    # tol = 1e-10
    
    # n = 2
    
    # htrap, errtrap = [], []
    # hsimp, errsimp = [], []
    
    # trapDone = False
    # simpDone = False
    
    # while not (trapDone and simpDone):
        
    #     h = (b-a)/n
        
    #     if not trapDone:
    #         I_trap = CompTrap(a,b,n,f)
    #         err_trap = abs(I_ex-I_trap)
    #         htrap.append(h)
    #         errtrap.append(err_trap)
    #         if err_trap/abs(I_ex) < tol:
    #             print('Composite trapezoidal rule converged with n =', n, 'h=', h)
    #             trapDone = True
                
    #     if not simpDone:
    #         I_simp = CompSimp(a,b,n,f)
    #         err_simp = abs(I_ex-I_simp)
    #         hsimp.append(h)
    #         errsimp.append(err_simp)
    #         if err_simp/abs(I_ex) < tol:
    #             print('Composite Simpson''s rule converged with n =', n, 'h=', h)
    #             simpDone = True
    #     n = n*2

    # slopeTrap = np.polyfit(np.log(htrap[-5:]),np.log(errtrap[-5:]),1)[0]
    # slopeSimp = np.polyfit(np.log(hsimp[-5:]),np.log(errsimp[-5:]),1)[0]
    # print('trapezoidal slope =', slopeTrap)
    # print('Simpson slope =', slopeSimp)
    
    # plt.loglog(htrap,errtrap,'ro--',label='trapezoidal')
    # plt.loglog(hsimp,errsimp,'bs--',label='Simpson')
    # plt.xlabel('h')
    # plt.ylabel('absolute error')
    # plt.legend()
    # plt.show()
    
    # nrefine = 5

# # number of panels: 1, 2, 4, ..., 2^nrefine
#     ntest = 2**np.arange(0,nrefine+1)

#     Rom = rombergTrap(a,b,nrefine,f)

# # column j is only filled from row j down
#     for j in range(0,nrefine+1):
#         print('absolute error in column', j, '= ', abs(I_ex-Rom[j:,j]))

#     err = abs(I_ex-np.diag(Rom))

#     print('absolute error of best approximation at each order = ', err)

#     plt.semilogy(ntest,err,'o-')
#     plt.xlabel('n')
#     plt.ylabel('absolute error')
#     plt.show()


def rombergTrap(a,b,nrefine,f):
# number of panels: 1, 2, 4, ..., 2^nrefine
    ntest = 2**np.arange(0,nrefine+1)

    nn = len(ntest)

    Rom = np.zeros((nn,nn))

# first column: composite trapezoidal rule with n panels
    for j in range(0,nn):
        Rom[j,0] = CompTrap(a,b,ntest[j],f)

# remaining columns: Richardson extrapolation
    for j in range(1,nn):
        for i in range(j,nn):
            Rom[i,j] = Rom[i,j-1] + (Rom[i,j-1]-Rom[i-1,j-1])/(4**j-1)

    return Rom


def CompTrap(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h

    I_trap = h*f(xnode[0])*1/2

    for j in range(1,n):
         I_trap = I_trap+h*f(xnode[j])
    I_trap= I_trap + 1/2*h*f(xnode[n])

    return I_trap

def CompSimp(a,b,n,f):
    h = (b-a)/n
    xnode = a+np.arange(0,n+1)*h

    I_simp = h*f(xnode[0])*1/3

    for j in range(1,n,2):
         I_simp = I_simp+h*f(xnode[j])*4/3
    for j in range(2,n-1,2):
         I_simp = I_simp+h*f(xnode[j])*2/3
    I_simp= I_simp + 1/3*h*f(xnode[n])

    return I_simp



driver()
