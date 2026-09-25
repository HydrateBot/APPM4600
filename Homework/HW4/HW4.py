"""
 This code implements Newton's method for finding the root of a
 scalar function.
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


# import libraries
import numpy as np
import matplotlib.pyplot as plt
        
def driver():
#f = lambda x: (x-2)**3
#fp = lambda x: 3*(x-2)**2
#p0 = 1.2

  f = lambda x: x**6 - x - 1
  fp = lambda x: 6*x**5 - 1
  p0 = 2
  p1 = 1
  alpha = 1.1347241384015194

  Nmax = 30
  tol = 1.e-14

  (p,pstar,info,it) = newton(f,fp,p0,tol, Nmax)
  
  print('the approximate root is', '%16.16e' % pstar)
  print('the error message reads:', '%d' % info)
  print('Number of iterations:', '%d' % it)
  
  newton_error = np.abs(p - alpha)
  
  print("\n Newton Method")
  print("k\t x_k \t\t\t\t error")
  
  for k in range(it + 1):
    print(k, '\t', "%16.16e" % p[k], '\t\t', "%16.16e" % newton_error[k])
  
  (xs, it, info) = secant(f,p0,p1,tol,Nmax)
  print('the approximate root is', '%16.16e' % xs[-1])
  print('the error message reads:', '%d' % info)
  print('Number of iterations:', '%d' % it)
  
  secant_error = np.abs(xs - alpha)

  print("\n Secant Method")
  print("k\t x_k \t\t\t\t error")

  for k in range(it + 1):
    print(k, '\t', "%16.16e" % xs[k], '\t\t', "%16.16e" % secant_error[k])

  newton_x_error = newton_error[:-1]
  newton_y_error = newton_error[1:]

  secant_x_error = secant_error[:-1]
  secant_y_error = secant_error[1:]
  
  newton_remove_zero = (newton_x_error > 0) & (newton_y_error > 0)
  secant_remove_zero = (secant_x_error > 0) & (secant_y_error > 0)
  
  newton_x_error = newton_x_error[newton_remove_zero]
  newton_y_error = newton_y_error[newton_remove_zero]
  secant_x_error = secant_x_error[secant_remove_zero]
  secant_y_error = secant_y_error[secant_remove_zero]

  plt.loglog(newton_x_error, newton_y_error, 'o-', label='Newton')
  plt.loglog(secant_x_error, secant_y_error, 's-', label='Secant')
  plt.xlabel('Error at iteration k')
  plt.ylabel('Error at iteration k+1')
  plt.legend()
  plt.show()
  
  newton_slope, newton_intercept = np.polyfit(np.log(newton_x_error), np.log(newton_y_error), 1)
  secant_slope, secant_intercept = np.polyfit(np.log(secant_x_error), np.log(secant_y_error), 1)
  
  print("Newton's method slope:", newton_slope)
  print("Secant method slope:", secant_slope)

def newton(f,fp,p0,tol,Nmax):
  """
  Newton iteration.
  
  Inputs:
    f,fp - function and derivative
    p0   - initial guess for root
    tol  - iteration stops when p_n,p_{n+1} are within tol
    Nmax - max number of iterations
  Returns:
    p     - an array of the iterates
    pstar - the last iterate
    info  - success message
          - 0 if we met tol
          - 1 if we hit Nmax iterations (fail)
     
  """
  p = np.zeros(Nmax+1);
  p[0] = p0
  for it in range(Nmax):
      p1 = p0-f(p0)/fp(p0)
      p[it+1] = p1
      if (abs(p1-p0) < tol):
          pstar = p1
          info = 0
          return [p,pstar,info,it]
      p0 = p1
  pstar = p1
  info = 1
  return [p,pstar,info,it]

def secant(f,p0,p1,tol,Nmax):
  xs = [p0, p1]
  
  for it in range(Nmax):
    
    x_prev = xs[-2]
    x_curr = xs[-1]
    
    fx_prev = f(x_prev)
    fx_curr = f(x_curr)
    
    x_next = x_curr - fx_curr * (x_curr - x_prev) / (fx_curr - fx_prev)
    xs.append(x_next)
    
    if abs(x_next - x_curr) < tol:
      info = 0
      return [np.array(xs), it + 1, info]
    
  info = 1
  return [np.array(xs), it, info]
        
driver()

