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
def f_example(x, y):
    return y - x

# or
# f = lambda x: y - x

def c_i(c0, k, t0, tn, h):
    """
    Equation 1 Dillenburger (2017)
    An exponential decay equation that calculates the cleavage of DMSP and
    formation of DMS at a certain time of the boil. 
    :param c0:  DMSP content at start of time span
    :param k:   Rate constant
    :param t0:  Start of time span
    :param tn:  End of time span
    :param h:   Time step
    :return:    DMSP content at end of time span
    """
    c = []
    for i in range(t0, tn+h, h):
        c.append(c0 * math.exp(-k * i))
    return c


def dc_i(t, c_i0=227, k=0.00025):
    """
    An alternative form of eq. 1 from Dillenburger (2017), presented as eq. 10 in Scheuren (2014).
    The right side of these two equations are equivelant. The difference is that the left side is now explicitly a
    differential equation that can be evaluated using the RK4 method.
    TODO make sure the parameters line up correctly with the paper.
    TODO implement this function in a way that it can be evaluated by the RK4 method.
    :param c_i: DMSP content
    :param t: Time
    :param k: Rate constant for chemical reaction.
    :return: DMSP content after a period of time in the boil.
    """
    return k * c_i0 * math.exp(-k * t)   # c_i is probably changing in each step, and it probably shouldn't for the eqn to be correct.


# Equation 4 Dillenburger (2017)
def f4(D_dot, L0, Ki, xi, ci, k, t):
    """
    Equation 4 from Dillenburger (2017) calculates 'the boiling of wort in a kettle through direct
    heating of tank walls in terms of the evaporation of DMS in water and teh simultaneous reproduction
    of DMSP'
    :param D_dot:   Steam flow (a constant or a variable?)
    :param L0:      Wort volume at point t0 (most likely a constant)
    :param Ki:      Volatility of component i (DMS) (most likely a constant)
    :param xi:      Dimethyl sulphide content
    :param ci:      Precursor content of component i (DMSP)
    :param k:       Rate constant 1/s
    :param t:       time
    :return:        The change in DMS content dependent on a process time
    """
    return



# RK-4 method
def rk4(f, t0, x0, h, n):
    """
    The parameters of this RK4 method have been altered to reflect those in the Scheuren/Dillenburger papers.
    :param f: The function to be analyzed with teh RK4 method.
    :param t0: Initial time for t in the form dx/dt. This corresponds to x in the form dy/dx.
    :param x0: Initial value of x in the form dx/dt. This corresponds to y in the form dy/dx.
    :param h: The time step for each iteration of the RK4 method.
    :param n: The total number of steps in the analysis.
    :return: A dataframe object containing the values of t and x at each step of the RK4 analysis for the duration of n.
    """
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
    # Inputs
    # TODO fix this so that it doesn't crash when NAN is entered
    print('Enter initial conditions:')
    t0 = float(input('t0 = '))
    x0 = float(input('x0 = '))
    print('Enter step size: ')
    h = float(input('h = '))
    print('Enter number of steps:')
    step = int(input('Number of steps = '))

    # RK4 method call dc_i
    data = rk4(dc_i, t0, x0, h, step)
    t = data["tn"]
    x = data["xn"]

    # RK4 method call example function
    # data = rk4(f_example, t0, x0, h, step)
    # t = data["tn"]
    # x = data["xn"]

    # plot the data
    plt.plot(t, x)
    plt.xlabel("t (in seconds)")
    plt.ylabel("x")
    plt.show()



