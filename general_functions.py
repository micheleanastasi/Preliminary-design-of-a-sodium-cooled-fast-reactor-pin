import numpy as np
import scipy as sp
from numpy import array,sqrt,pi

from material_properties import *
from design_specifications import *

#### GENERAL FUNCTIONS ####

def equation_temp_solver(sy_funct, guess):
    """
    Give result from interaction according to one variable: temperature
    :param sy_funct: function to be iterated in Sympy format and temp as variable
    :param guess: initial value
    :return:
    """
    res = (sy_funct)**2
    res_l = sy.lambdify(temp,res)
    out = sp.optimize.minimize((res_l),guess,tol=10e-3)
    return out.x[0]



#### *********************************** POWER FUNCTIONS ********************************* ####

## LINEAR POWER DISTRIBUTION FUNCTIONS - INTERPOLATED ##
def interpolated_power(z):
    unit = pin_top_pos / 10  #0.085
    peak_factor = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * power_lin_max

    discrete_domain = np.arange(pin_bottom_pos + unit / 2, pin_top_pos - unit / 2, unit)
    yy_power = np.zeros_like(discrete_domain)
    for i in range(0, len(discrete_domain)):
        value = int(discrete_domain[i] / unit)
        yy_power[i] = discrete_power[value]

    coeff = np.polyfit(discrete_domain, yy_power, 9)
    poly = np.poly1d(coeff)
    y_new = poly(z)
    return y_new

def integral_cosine_shape_power(z):
    area = sp.integrate.quad(interpolated_power, pin_bottom_pos, z)
    return area[0]


## LINEAR POWER DISTRIBUTION FUNCTIONS - STEP FUNCTIONS ##
#
## approx or not: change True or False!!
def power_lin_distribution(z,approx=True):
    unit = pin_top_pos / 10  #0.085
    peak_factor = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * power_lin_max
    value = int(z/unit)
    output = discrete_power[value]

    if approx == True:
        return interpolated_power(z)
    else:
        return output

def integral_power_lin_distr(z,approx=True):
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

    if approx == True:
        return integral_cosine_shape_power(z)
    else:
        return area