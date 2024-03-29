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
import pandas as pd
from matplotlib import pyplot as plt


###### Function to be solved ######
# Example problem
def f_example(x, y):
    """
    A simple example problem to confirm that the RK4 method is working correctly.
    """
    return y - x


def eqn_1(c_i_0, k, t, step):
    """
    Equation 1 from Scheuren (2014) Thermodynamic Validation of Wort Boiling Systems.

    Calculation of the cleavage of DMSP and formation of DMS
    Paper describes it as a differential equation but it does not resemble one:
    c_i,1 = c_i,0 * e**(-k * t_1)

    NOTE 1: Scheuren uses 2.71 for e in the Excel file he gave me,
        This returns slightly different values than when using math.exp()

    NOTE 2: The 60 is to convert an input time of minutes to seconds.
        This does not appear in the paper, but is included in the Excel document from Scheuren. 

    :param c_i_0: Initial DMSP content.
    :param k: Rate constant.
    :param t: Given instant in time.
    :param step: Time interval for each calculation.
    :return: List with DMSP concentration for each time step in the given range. 
    """
    c_i = []
    for i in range(0, t + step, step):
        c_i.append(c_i_0 * math.exp(-k * 60 * i))
        # c_i.append(c_i_0 * 2.71**(-k * 60 * i))   # NOTE: Using 2.71 instead of math.exp will cause unittest to fail.
    
    return c_i
    

D = [
    0,
    0.005,
    0.01,
    0.015,
    0.02,
    0.025,
    0.03,
    0.035,
    0.04,
    0.045,
    0.05,
    0.055,
    0.06,
    0.065,
    0.07,
    0.075,
]

def eqn_4():
    """
    Equation 4 from Scheuren (2014) Thermodynamic Validation of Wort Boiling Systems.

    Equation models the evaporation of DMS from water kettle boiling without consideration of DMSP to DMS conversion.

    (dx_i)/(dt) = (-D_dot * x_i * (K_i - 1)) / L

    :param x_i: DMS content
    :param D_dot: Overall evaporation
    :param K_i: Volatility or fugacity
    :param L: Volume
    :return: DMS content????
    """
    # result = [
    #     500,
    #     343.3215466,
    #     235.2933208,
    #     160.9486327,
    #     109.8817788,
    #     74.87143927,
    #     50.91550719,
    #     34.55560533,
    #     23.40520678,
    #     15.82059719,
    #     10.67186692,
    #     7.183823317,
    #     4.825685608,
    #     3.234750732,
    #     2.16367061,
    #     1.444108951
    # ]
    result = 1
    return result

def eqn_13(t, x_i):
    """
    Equation 13 from Scheuren (2014) Thermodynamic Validation of Wort Boiling Systems.

    Equation 4 from Dillenburger (2017) calculates 'the boiling of wort in a kettle through direct
    heating of tank walls in terms of the evaporation of DMS in water and teh simultaneous reproduction
    of DMSP'
    
    Local variables go against python convention to keep them consistent with the variables used in the papers.
    :param x_i: Dimethyl sulphide content
    :param t: time
    :return: The change in DMS content dependent on a process time

    NOTE: If using Dillenburger values, make sure to alter the x0 input for the RK4 Method (see belowe)
    Notes on results:
    - Using all Dillenburger constants and x0 = 500 yields a boil time of 30 minutes (1760ish seconds) to reach 100 ug/l of DMS
    - Using all Dillenburger constants and x0 = 273 yields a boil time of 20 minutes (1209ish seconds) to reach 100 ug/l of DMS
    - Using constants from Scheuren's Excel where available and x0 = 500 yields a boil time of 43 minutes (2590 seconds)
    """
    # Steam flow in L/s
    D_dot = 1.1     # From Dillenburger (2017)

    # Initial wort volume (boil start)
    L_0 = 73200     # From Dillenburger (2017)

    # Volatility
    K_i = 78        # From Dillenburger (2017) (was 76 before refactor, not sure where 76 came from)

    # Rate constant
    k = 0.00130809768004509     # From Scheuren's Excel
    # k = 0.00025     # From Dillenburger (2017)
    
    # Dimethyl sulphide precursor (DMSP) in micrograms/L
    c_i0 = 500      # From Scheuren's Excel
    # c_i0 = 227      # From Dillenburger (2017)

    return -(D_dot/L_0) * (K_i * x_i - x_i - c_i0) + k * c_i0 * math.exp(-k * t)


def rk4(f, t0, x0, h, n):
    """
    Runge-Kutta 4th Order method for solving a differential equation. 

    :param f: The function to be analyzed with teh RK4 method.
    :param t0: Initial time for t in the form dx/dt. This corresponds to x in the form dy/dx.
    :param x0: Initial value of x in the form dx/dt. This corresponds to y in the form dy/dx.
    :param h: The time step for each iteration of the RK4 method.
    :param n: The total number of steps in the analysis.
    :return: A dictionary with tn as each key and xn as the corresponding value.
    :return: A list with the header names.
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
        data[round(tn, 2)] = xn
        xn = xn + k

    return data, headers


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

    # NOTE: When modeling eqn 4 from Dillenburger (2017), input the following:
    # NOTE: The following inputs came from an interpretation of the 2017 Dillenburger paper.
    # t0 = 0  # (Start time)
    # x0 = 273    # (DMS content)
    # h = 1   # (Step size of 1 second)
    # steps = 3600    # (3600 seconds is 1 hr)

    ##### RK4 method call f4 ######
    data, headers = rk4(eqn_13, t0, x0, h, steps)
    df = pd.DataFrame(data.items(), columns=headers)
    print("\nRK4 DataFrame: \n", df)
    t = df["tn"]
    x = df["xn"]

    ###### plot the data #####
    plt.axhline(y=100, color='red')
    plt.plot(t, x)
    plt.xlabel("t (in seconds)")
    plt.ylabel("DMS Content (x) in ug/l")
    plt.show()

    ###### RK4 method call example function ######
    ###### This is used to test/confirm that the RK4 method works correctly ######
    # t0 = 0
    # x0 = 2
    # h = 0.1
    # steps = 20
    # data, headers = rk4(f_example, t0, x0, h, steps)
    # print("data:\n", data)
    # print("headers:\n", headers)
    # df = pd.DataFrame(data.items(), columns=headers)
    # print("\nRK4 DataFrame: \n", df)
    # t = df["tn"]
    # x = df["xn"]

    # plt.plot(t, x)
    # plt.xlabel("t (in seconds)")
    # plt.ylabel("DMS Content (x) in ug/l")
    # plt.show()

