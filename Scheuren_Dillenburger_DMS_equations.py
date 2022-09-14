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

    return


# Equation 4 Dillenburger (2017)
def f(x, y):

    return



# RK-4 method
def rk4(x0, y0, h, n):
    # Calculating step size
    # h = (xn - x0) / n

    # Calculate initial conditions
    # xn = n * h
    yn = y0
    data = {x0: yn}
    n = n + 1
    headers = ["xn", "yn"]

    for i in range(n):
        xn = i * h
        k1 = h * (f(xn, yn))
        k2 = h * (f((xn + h / 2), (yn + k1 / 2)))
        k3 = h * (f((xn + h / 2), (yn + k2 / 2)))
        k4 = h * (f((xn + h), (yn + k3)))
        k = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        data[xn] = yn
        yn = yn + k

    df = pd.DataFrame(data.items(), columns=headers)
    print(df)
    return df


if __name__ == "__main__":
    # Inputs
    print('Enter initial conditions:')
    x0 = float(input('x0 = '))
    y0 = float(input('y0 = '))

    print('Enter step size: ')
    h = float(input('h = '))

    print('Enter number of steps:')
    step = int(input('Number of steps = '))

    # RK4 method call
    data = rk4(x0, y0, h, step)
    # print(data.to_numpy()[0][0])
    # print(data["xn"])
    # print(data["yn"])

    x = data["xn"]
    y = data["yn"]

    # matplotlib example
    # x = np.linspace(0, 2 * np.pi, 200)
    # y = np.sin(x)

    fig, ax = plt.subplots()
    ax.plot(x, y)
    plt.show()



