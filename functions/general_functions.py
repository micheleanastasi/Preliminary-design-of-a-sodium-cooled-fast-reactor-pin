"""
General function involving functions useful for more stuff
Power distribution function and its integral
"""
import numpy as np
import scipy as sp
from numpy import array,pi

from functions.general_properties import *


#### ************************************************************************************************************** ####
#### ********************************************** GENERAL FUNCTIONS ********************************************* ####
#### ************************************************************************************************************** ####

def sy_equation_solver(sy_funct, guess):
    """
    Give result from interaction according to one variable: temperature
    :param sy_funct: function to be iterated in Sympy format and temp as variable
    :param guess: initial value
    :return:
    """
    res = sy_funct ** 2
    res_l = sy.lambdify(temp,res) # conversione sympy --> def python to be used with minimize (scypy)
    out = sp.optimize.minimize(res_l, guess, tol=10e-3)
    return out.x[0] # returning solved equation (temperature)

def fun_equation_solver(funct, guess):
    res = lambda t : funct(t) ** 2
    out = sp.optimize.minimize(res, guess, tol=10e-3)
    return out.x[0]

def volume_calc(d,h):
    return (h*pi*d**2)/4


def burnup_calc(days):
    vol = 0.25*3.1415*pin_column_height*fuel_d_outer**2
    mass = fuel_density * vol
    factor = power_lin_max*pin_column_height/mass
    return factor*days



#### ************************************************************************************************************** ####
#### *********************************************** POWER FUNCTIONS ********************************************** ####
#### ************************************************************************************************************** ####

## INTERPOLATED - LINEAR POWER DISTRIBUTION FUNCTIONS ##
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

## INTERPOLATED - LINEAR POWER INTEGRAL (up to z) ##
def integral_interpolated_power(z):
    area = sp.integrate.quad(interpolated_power, pin_bottom_pos, z)
    return area[0]



## STEP FUNCTION - LINEAR POWER DISTRIBUTION ##
def step_power(z):
    unit = pin_top_pos / 10  #0.085
    peak_factor = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * power_lin_max
    value = int(z/unit)
    output = discrete_power[value]
    return output

## STEP FUNCTION - LINEAR POWER INTEGRAL (up to z) ##
def integral_step_power(z):
    discrete_power = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498]) * power_lin_max
    area = 0
    unit = pin_top_pos / 10  # 0.085
    value = int(z / unit)
    total = np.sum(discrete_power) * unit

    if value == 10:  # if z = 0.85  then get all the power (pin)
        area = total
    else:
        for i in range(0, value + 1):
            if i == value:  # if we are in the last node, then not giving whole of the power but only a piece of the step
                area += discrete_power[i] * (z - i * unit)
            else:  # otherwise give the area (power) of the steps
                area += discrete_power[i] * unit
    return area


## POWER FUNCTIONS CHOICE ##
def power_lin_distribution(z,approx=True):
    """
    :param z: in m
    :param approx: interpolated (true), step function (false)
    :return: Return linear power distribution
    """
    if approx:
        return interpolated_power(z)
    else:
        return step_power(z)


def integral_power_lin_distr(z,approx=True):
    """
    Calculate integral, from bottom pos of pin (0) to z, of the power distribution function above!
    :return:  power of the pin up to z, in [W]
    """
    if approx:
        return integral_interpolated_power(z)
    else:
        return integral_step_power(z)



#### ************************************************* MECHANICS ************************************************** ####

def contact_pressure(fuel_r_out, clad_r_in, clad_r_out, fuel_poisson_ratio, fuel_young_modulus, clad_poisson_ratio, clad_young_modulus, fuel_r_void=0):
    interference = fuel_r_out - clad_r_in
    if interference >= 0:
        p = (interference / clad_r_in) / (((clad_r_in ** 2 + clad_r_out ** 2) / (clad_r_out ** 2 - clad_r_in ** 2) + 1 / clad_poisson_ratio) / clad_young_modulus + ((clad_r_in ** 2 + fuel_r_void ** 2) / (clad_r_in ** 2 - fuel_r_void ** 2) - 1 / fuel_poisson_ratio) / fuel_young_modulus)
    else:
        p = 0
    return p


