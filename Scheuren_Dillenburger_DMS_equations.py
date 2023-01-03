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
from cmath import exp
import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


###### Function to be solved ######
# Example problem
def f_example(x, y):
    """
    A simple example problem to confirm that the RK4 method is working correctly.
    """
    return y - x


# Equation 1 from Scheuren (2014) Thermodynamic Validation of Wort Boiling Systems
def eqn_1(c_i_0, k, t, step):
    """
    Calculation of the cleavage of DMSP and formation of DMS
    Paper describes it as a differential equation but it does not resemble one:
    c_i,1 = c_i,0 * e**(-k * t_1)
    Note 1: Scheuren uses 2.71 for e in the Excel file he gave me,
        This returns slightly different values than when using math.exp()
    Note 2: The 60 is to convert an input time of minutes to seconds
    """
    c_i = []
    for i in range(0, t + step, step):
        c_i.append(c_i_0 * math.exp(-k * 60 * i))
    
    return c_i
    

# Equation 4 Dillenburger (2017)
def eqn_4(t, x_i):
    """
    Equation 4 from Dillenburger (2017) calculates 'the boiling of wort in a kettle through direct
    heating of tank walls in terms of the evaporation of DMS in water and teh simultaneous reproduction
    of DMSP'
    Local variables go against python convention to keep them consistent with the equation in the paper.
    :param x_i: Dimethyl sulphide content
    :param t:   time
    :return:    The change in DMS content dependent on a process time
    """
    D_dot = 1.1     # Steam flow in L/s
    L_0 = 73200     # Initial wort volume (boil start)
    K_i = 76        # Volatility
    k = 0.00130809768004509     # Rate constant
    c_i0 = 500      # Dimethyl sulphide precursor (DMSP) in micrograms/L

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
    return df


if __name__ == "__main__":


    # For eqn 1
    c_i_0 = 500    # From Scheuren's Excel 
    k = 0.00130809768004509     # From Scheuren's Excel
    c_i = eqn_1(c_i_0, k, 60, 5)
    print('c_i:\n', c_i)

    # For RK4
    # Inputs
    t0 = 0
    x0 = 500
    h = 1
    steps = 3600

    # Note: When modeling eqn 4, input the following:
    # t0 = 0 (Start time)
    # x0 = 273 (DMS content)
    # h = 1 (Step size of 1 second)
    # Number of Steps = 3600 (3600 seconds is 1 hr)

    ##### RK4 method call f4 ######
    data = rk4(eqn_4, t0, x0, h, steps)
    print("\nRK4 DataFrame: \n", data)
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
    plt.ylabel("DMS Content (x) in ug/l")
    plt.show()
