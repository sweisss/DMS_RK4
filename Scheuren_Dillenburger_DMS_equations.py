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
import math
import numpy as np
import pandas as pd
from tabulate import tabulate
from matplotlib import pyplot as plt


# RK-4 method python program
# from codesansar.com

# function to be solved
# Example problem
# def f(x, y):
#     return y - x

# or
# f = lambda x: y - x

def c_i(c0, k, t0, tn, h):
    """
    Equation 1 Dillenburger (2017)
    An exponential decay equation that calculates the DMSP content
    at a certain time of the evaporation process.
    :param c0:  DMSP content at start of time span
    :param k:   Rate constant
    :param t0:  Start of time span
    :param tn:  End of time span
    :param h:   Time step
    :return:    DMSP content at end of time span
    """
    c = []
    for i in range(t0, tn+h, h):
        c.append(c0 * math.exp(-k*i))
    return c


# Equation 4 Dillenburger (2017)
def f(t, x):
    return x - t



# RK-4 method
def rk4(t0, x0, h, n):
    # Calculating step size
    # h = (xn - x0) / n

    # Calculate initial conditions
    # xn = n * h
    xn = x0
    data = {t0: xn}
    n = n + 1
    headers = ["tn", "xn"]

    for i in range(n):
        tn = i * h
        k1 = h * (f(tn, xn))
        k2 = h * (f((tn + h / 2), (xn + k1 / 2)))
        k3 = h * (f((tn + h / 2), (xn + k2 / 2)))
        k4 = h * (f((tn + h), (xn + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        data[tn] = xn
        xn = xn + k

    df = pd.DataFrame(data.items(), columns=headers)
    print(df)
    return df


if __name__ == "__main__":
    # Inputs -- fix this so that it doesn't crash when NAN is entered
    print('Enter initial conditions:')
    t0 = float(input('t0 = '))
    x0 = float(input('x0 = '))
    print('Enter step size: ')
    h = float(input('h = '))
    print('Enter number of steps:')
    step = int(input('Number of steps = '))

    # RK4 method call
    data = rk4(t0, x0, h, step)
    t = data["tn"]
    x = data["xn"]

    # plot the data
    fig, ax = plt.subplots()
    ax.plot(t, x)
    plt.show()



