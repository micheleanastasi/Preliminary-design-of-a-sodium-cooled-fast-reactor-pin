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
    Give result from interaction according to one variable: temperature
    :param sy_funct: function to be iterated in Sympy format and temp as variable
    :param guess: initial value
    :return:
    """
    res = (sy_funct)**2
    res_l = sy.lambdify(temp,res)
    out = sp.optimize.minimize((res_l),guess,tol=10e-3)

    #print(out)
    return out.x[0]

def plotting(funct,name=None,legend=None,x_label=None,y_label=None,x_lim=None,y_lim=None,plot=True):
    """
    :param funct: function to be plotted
    :return: plot
    """
    xx = np.linspace(pin_bottom_pos, pin_top_pos, 100)
    yy = np.zeros(100)

    for i in range(0,len(xx)):
        yy[i] = funct(xx[i])
    plt.figure(name)
    if legend != None:
        plt.plot(xx,yy,label = legend)
    else:
        plt.plot(xx, yy)
    if x_label is not None:
        plt.xlabel(x_label)
    if y_label is not None:
        plt.ylabel(y_label)
    if x_lim is not None:
        plt.xlim(x_lim)
    if y_lim is not None:
        plt.ylim(y_lim)
    if name is not None:
        plt.title(name)
    if legend is not None:
        plt.legend()
    plt.grid()
    if plot is True:
        plt.show()
    return None

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

#### END OF POWER FUNCT ####

