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
from matplotlib import pyplot as plt


###### Function to be solved ######
# Example problem
def f_example(x, y):
    return y - x


# Equation 4 Dillenburger (2017)
def f4(t, x_i):
    """
    Equation 4 from Dillenburger (2017) calculates 'the boiling of wort in a kettle through direct
    heating of tank walls in terms of the evaporation of DMS in water and teh simultaneous reproduction
    of DMSP'
    Local variables go against python convention to keep them consistent with the equation in the paper.
    :param x_i:      Dimethyl sulphide content
    :param t:       time
    :return:        The change in DMS content dependent on a process time
    """
    D_dot = 1.1     # Steam flow in L/s
    L_0 = 73200     # Initial wort volume (boil start)
    K_i = 78        # Volatility
    k = 0.00025     # Rate constant
    c_i0 = 227      # Dimethyl sulphide precursor (DMSP) in micrograms/L

    return -(D_dot/L_0) * (K_i * x_i - x_i - c_i0) + k * c_i0 * math.exp(-k * t)


###### RK-4 method ######
def rk4(f, t0, x0, h, n):
    """
    The parameters of this RK4 method have been altered to reflect those in the Scheuren/Dillenburger papers.
    :param f: The function to be analyzed with teh RK4 method.
    :param t0: Initial time for t in the form dx/dt. This corresponds to x in the form dy/dx.
    :param x0: Initial value of x in the form dx/dt. This corresponds to y in the form dy/dx.
    :param h: The time step for each iteration of the RK4 method.
    :param n: The total number of steps in the analysis.
    :return: A dataframe object containing the values of t and x at each step of the RK4 analysis for the duration of n.
    TODO pull the df part out of this and just return the dictionary. The df part should be a separate function.
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
    t0 = 0
    x0 = 273
    h = 1
    steps = 3600

    # Note: When modeling eqn 4, input the following:
    # t0 = 0 (Start time)
    # x0 = 273 (DMS content)
    # h = 1 (Step size of 1 second)
    # Number of Steps = 3600 (3600 seconds is 1 hr)

    ###### RK4 method call f4 ######
    data = rk4(f4, t0, x0, h, steps)
    t = data["tn"]
    x = data["xn"]

    ###### RK4 method call example function ######
    # data = rk4(f_example, t0, x0, h, steps)
    # t = data["tn"]
    # x = data["xn"]

    ###### plot the data #####
    plt.axhline(y=100, color='red')
    plt.plot(t, x)
    plt.xlabel("t (in seconds)")
    plt.ylabel("x")
    plt.show()
