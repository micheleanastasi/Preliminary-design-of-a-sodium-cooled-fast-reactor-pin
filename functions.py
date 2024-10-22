import numpy as np
import scipy as sp
from numpy import array
import matplotlib.pyplot as plt
import sympy as sy

from material_properties import *
from design_specifications import *



#### POWER FUNCT ####

## Symbols used ##
_z = sy.Symbol('z')

## LINEAR POWER DISTRIBUTION FUNCTION ##
def power_lin_distribution(z):
    unit = pin_top_pos / 10  #0.085
    peak_factor = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * power_lin_max
    #discrete_pos = np.arange(_unit / 2, pin_top_pos - _unit / 2, _unit)
    value = int(z/unit)
    print(value)

    output = discrete_power[value]
    return output

def integral_power_lin_distr():
    """
    :return:  power of the pin in W
    """
    sum = 0
    unit = pin_top_pos / 10  #0.085
    discrete_power = array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498])*power_lin_max
    for i in range(0,len(discrete_power)-1):
        sum += discrete_power[i] * unit

    return sum

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
    plt.plot(xx,yy)
    plt.xlabel("Fuel pin axis in [m]")
    plt.ylabel("Linear power in [W/m]")
    plt.title("Linear power distribution")
    plt.show()
    return None

#### TEMPERATURE FUNCTIONS ####

## TEMP PROFILE OF COOLANT ALONG Z AXIS ##
def temp_coolant(z):
    """
    :param z: position choice
    :return: Temp of coolant at point T
    """
    # energy balance: int(Cp, Tin, T) = int(q', bottom_pin, z) / mass_flow ... Tout our goal.




    return None

temp_coolant(0.5)