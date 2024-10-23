import numpy as np
import scipy as sp
from numpy import array,sqrt,pi
import matplotlib.pyplot as plt
import sympy as sy

from material_properties import *
from design_specifications import *

#### GENERAL FUNCTIONS ####

def iterative_solver(sy_funct, guess):
    """
    Give result from interaction according to one variable
    :param sy_funct: function to be iterated in Sympy format
    :param guess: initial value
    :return:
    """
    res = (sy_funct)**2
    res_l = sy.lambdify(temp,res)
    out = sp.optimize.minimize((res_l),guess,tol=10e-3)

    return out.x[0]

#### END OF GENERAL FUNCTIONS ####

#### POWER FUNCT ####

## LINEAR POWER DISTRIBUTION FUNCTION ##
def power_lin_distribution(z):

    unit = pin_top_pos / 10  #0.085
    peak_factor = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * power_lin_max
    #discrete_pos = np.arange(_unit / 2, pin_top_pos - _unit / 2, _unit)
    value = int(z/unit)

    output = discrete_power[value]

    return output


def integral_power_lin_distr(z):
    """
    Calculate integral, from bottom pos of pin (0) to z, of the power distribution function above!
    :return:  power of the pin up to z, in [W]
    """
    discrete_power = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498])*power_lin_max
    
    area = 0
    unit = pin_top_pos / 10  #0.085
    value = int(z/unit)
    
    total = np.sum(discrete_power)*unit
    
    if value == 10:     # if z = 0.85  then get all the power (pin)
        area = total
    else:
        for i in range(0,value+1):
            if i == value:  # if we are in the last node, then not giving whole of the power but only a piece of the step
                area += discrete_power[i] * (z - i*unit)
            else: # otherwise give the area (power) of the steps
                area += discrete_power[i]*unit
    return area


## Plot function for power distr.
def plot_power_distribution():
    """
    Simply print the plot of power distribution, only fuel pin dimensions (0 to 0.85 m)

    :param res: function to be plotted (from bottom to top of the pin)
    :return: plot
    """
    num_range = 100
    xx = np.linspace(pin_bottom_pos, pin_top_pos, num_range)
    yy = np.zeros(num_range) # initialized
    for i in range(0,len(yy)):
        yy[i] = power_lin_distribution(xx[i])
    plt.figure("Linear power distribution")
    plt.plot(xx,yy)
    plt.xlabel("Fuel pin axis in [m]")
    plt.ylabel("Linear power in [W/m]")
    plt.title("Linear power distribution")
    plt.ylim([0,40000])
    plt.grid()
    plt.show()
    return None

def plot_integral_power_distribution():
    """
    Simply print the plot of power distribution, only fuel pin dimensions (0 to 0.85 m)

    :param res: function to be plotted (from bottom to top of the pin)
    :return: plot
    """
    num_range = 100
    xx = np.linspace(pin_bottom_pos, pin_top_pos, num_range)
    yy = np.zeros(num_range) # initialized
    for i in range(0,len(yy)):
        yy[i] = integral_power_lin_distr(xx[i])
    plt.figure("Linear power distribution")
    plt.plot(xx,yy)
    plt.xlabel("Fuel pin axis in [m]")
    plt.ylabel("Linear power in [W/m]")
    plt.title("Linear power distribution")
    plt.ylim([0,30000])
    plt.grid()
    plt.show()
    return None

#### END OF POWER FUNCT ####


plot_power_distribution()
plot_integral_power_distribution()