import numpy as np

def f(x, y, z):
    return x**2 + 4*y**2 + 4*z**2 - 16

def grad_f(x, y, z):
    fx = 2*x
    fy = 8*y
    fz = 8*z

    return np.array([fx, fy, fz])

x = 1.0
y = 1.0
z = 1.0

tol = 1e-10
max_iter = 20

points = []

print("Iteration results:")
print(" n          x_n            y_n            z_n            f(x,y,z)          d")
print("-" * 95)

for n in range(max_iter):
    
    points.append(np.array([x, y, z]))
    
    fval = f(x, y, z)
    grad = grad_f(x, y, z)

    fx = grad[0]
    fy = grad[1]
    fz = grad[2]

    d = fval / (fx**2 + fy**2 + fz**2)

    print(f"{n:2d}   "
          f"{x:14.10f} "
          f"{y:14.10f} "
          f"{z:14.10f} "
          f"{fval:16.10e} "
          f"{d:16.10e}")

    if abs(fval) < tol:
        break

    x_new = x - d*fx
    y_new = y - d*fy
    z_new = z - d*fz

    x = x_new
    y = y_new
    z = z_new

solution = np.array([x, y, z])

print("\nApproximate solution:")
print(f"x = {x:.12f}")
print(f"y = {y:.12f}")
print(f"z = {z:.12f}")

print("\nCheck:")
print(f"x^2 + 4y^2 + 4z^2 = {x**2 + 4*y**2 + 4*z**2:.12f}")

print("\nQuadratic convergence test:")
print(" n           error e_n               e_(n+1) / e_n^2")
print("-" * 65)

errors = []

for point in points:
    
    error = np.linalg.norm(point - solution)
    
    errors.append(error)
    
for n in range(len(errors) - 1):
    
    e_n = errors[n]
    e_next = errors[n + 1]
    
    if e_n > 0:
        ratio = e_next / e_n**2
        print(f"{n:2d}       "      f"{e_n:16.10e}        "         f"{ratio:16.10e}")