"""
Written by: Seth Weiss
Started on: 9/12/22
Description: A program that uses the Runge Kutta method (RK4) of numerical analysis
to solve the differential equations presented in papers by H. Scheuren and M. Dillenburger.
This program starts by lifting the RK4 example from Geeks for Geeks written by Arpit Agarwal:
https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
This program also used the code from the following link as a starting point:
https://www.codesansar.com/numerical-methods/runge-kutta-fourth-order-rk4-python-program.htm
"""
import numpy as np
import pandas as pd
from tabulate import tabulate
from matplotlib import pyplot as plt

# # Python program to implement Runge Kutta method
# # From G4G
# # A sample differential equation "dy / dx = (x - y)/2"
# def dydx(x, y):
#     return ((x - y) / 2)
#
#
# # Finds value of y for a given x using step size h
# # and initial value y0 at x0.
# def rungeKutta(x0, y0, x, h):
#     # Count number of iterations using step size or
#     # step height h
#     n = (int)((x - x0) / h)
#     # Iterate for number of iterations
#     y = y0
#     for i in range(1, n + 1):
#         "Apply Runge Kutta Formulas to find next value of y"
#         k1 = h * dydx(x0, y)
#         k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
#         k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
#         k4 = h * dydx(x0 + h, y + k3)
#
#         # Update next value of y
#         y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
#
#         # Update next value of x
#         x0 = x0 + h
#     return y


# RK-4 method python program
# from codesansar.com

# function to be solved
def f(x, y):
    return y - x


# or
# f = lambda x: x+y

# RK-4 method
def rk4(x0, y0, h, n):
    # Calculating step size
    # h = (xn - x0) / n

    # Calculate initial conditions
    # xn = n * h
    yn = y0
    data = {x0: yn}
    n = n+1
    headers = ["xn", "yn"]

    for i in range(n):
        xn = i * h
        k1 = h * (f(xn, yn))
        k2 = h * (f((xn + h / 2), (yn + k1 / 2)))
        k3 = h * (f((xn + h / 2), (yn + k2 / 2)))
        k4 = h * (f((xn + h), (yn + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        # print('%.4f\t%.4f' % (xn, yn))
        # print('-------------------------')
        data[xn] = yn
        yn = yn + k

    # print('\nAt x=%.4f, y=%.4f' % (xn, yn))
    print(pd.DataFrame(data.items(), columns=headers))

if __name__ == "__main__":

    # # This code is contributed by Prateek Bhindwar
    # # Driver method
    # x0 = 0
    # y = 1
    # x = 2
    # h = 0.2
    # print('The value of y at x is:', rungeKutta(x0, y, x, h))



    # Inputs
    print('Enter initial conditions:')
    x0 = float(input('x0 = '))
    y0 = float(input('y0 = '))

    print('Enter step size: ')
    h = float(input('h = '))

    print('Enter number of steps:')
    step = int(input('Number of steps = '))

    # RK4 method call
    rk4(x0, y0, h, step)

