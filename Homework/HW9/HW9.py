import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

def driver():
    f = lambda x: np.sin(x)
    
    p33 = lambda x: (x - (7/60)*x**3) / (1 + x**2/20)
    p42 = lambda x: (x - (7/60)*x**3) / (1 + x**2/20)
    p24 = lambda x: (x / (1 + x**2/6 + (7/360)*x**4))
    
    t = lambda x: x - x**3/6 + x**5/120
    
    x = np.linspace(0, 5, 1000)
    
    error33 = np.abs(f(x) - p33(x))
    error42 = np.abs(f(x) - p42(x))
    error24 = np.abs(f(x) - p24(x))
    errorTaylor = np.abs(f(x) - t(x))

    # Plot the errors
    plt.figure(figsize=(10, 6))

    plt.plot(x, error33, label="Pade (3,3)")
    plt.plot(x, error42, "--", label="Pade (4,2)")
    plt.plot(x, error24, label="Pade (2,4)")
    plt.plot(x, errorTaylor, label="Maclaurin Degree 6")

    plt.xlabel("x")
    plt.ylabel("Absolute Error")
    plt.title("Pade Approximation Errors for sin(x)")

    plt.legend()
    plt.grid(True)
    plt.show()

    # Print maximum errors
    print("Maximum Error Pade (3,3):", np.max(error33))
    print("Maximum Error Pade (4,2):", np.max(error42))
    print("Maximum Error Pade (2,4):", np.max(error24))
    print("Maximum Error Maclaurin:", np.max(errorTaylor))


driver()
