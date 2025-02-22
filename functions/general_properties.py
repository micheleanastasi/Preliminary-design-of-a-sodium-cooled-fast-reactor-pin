"""
Material properties deriving from homework.pdf data
"""

import sympy as sy
from sympy import exp
import numpy as np


#### DATA GUESS ####
clad_thickness_0 = 0.525e-3 # m
extra_pin_len = 0.70    # m - little diameter expansion then (whereas length exp neglected!) (HP CONS) !

#### ************************************************************************************************************** ####
#### ************************************************ DESIGN SPECS ************************************************ ####
#### ************************************************************************************************************** ####


#### ************** POWER DISTRIBUTION ************* ####
power_lin_max = 38700 # W/m - linear power @ peak factor node (@ 0.3825 m over the bottom)
peak_factor = np.array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498])
power_lin_avg=np.average(power_lin_max*peak_factor)
neu_max_flux = 6.1e15 # n/cm^2/sec


#### ************ THERMO HYDRAULICS ************* ####
cool_temp_inlet = 395 + 273.15 # K
cool_mass_flow = 0.049 # kg/s
cool_press_inlet = 0.1e6 # Pa
f_pitch = 0.5 # volume control is a triangle
pin_pitch = 8.275e-3 # m


#### ***************** FUEL PIN SPECIFICATIONS ******************* ####
temp_in = 293.15 # K

pin_bottom_pos = 0  # m
pin_top_pos = .850  # m
pin_column_height = 0.85 # m

clad_d_outer = 6.55e-3 # m
clad_d_inner = clad_d_outer - clad_thickness_0*2

fuel_d_outer = 5.42e-3 #m
fuel_height = 7e-3 # m

initial_delta_gap = clad_d_outer/2 - clad_thickness_0 - fuel_d_outer/2

fill_gas_temp_in = 293.15 # K - initial temperature
fill_gas_press_in = 1e5 # Pa

#### ************** POWER DISTRIBUTION ************* ####
power_lin_max = 38700 # W/m - linear power @ peak factor node (@ 0.3825 m over the bottom)

neu_max_flux = 6.1e15 # n/cm^2/sec




#### ************************************************************************************************************** ####
#### *********************************************** MATERIAL SPECS *********************************************** ####
#### ************************************************************************************************************** ####

## Symbols used from sympy - conversioni temperatura
temp = sy.Symbol('T[K]')
temp_f = sy.Symbol('T[°F]')
Pe_cool = sy.Symbol('Pe')
x_om = sy.Symbol('x[O/M]')
pu_conc = sy.Symbol('[Pu]')
por = sy.Symbol('P')
b_up = sy.Symbol('B-up')


#### ****************** COOLANT PROPERTIES ********************* ####

cool_T_melting_atm = 98 + 273.15 # K @ p = 0.1 Mpa
cool_T_boiling_atm = 882 + 273.15 # K @ p = 0.1 Mpa

## Coolant specific heat - input in K - output in J/kg/K
cool_spec_heat = 1608 - 0.7481*temp + 3.929e-4*temp ** 2

## Coolant density - input in °F!!!!!! - output in kg/m^3
cool_density__fahr = 954.1579 + temp_f*( temp_f*(temp_f*0.9667e-9 - 0.46e-5) - 0.1273534 )

## Coolant dyn viscosity - input in K - output in Pa*s
cool_dyn_viscosity = exp( 813.9/temp - 2.530 ) * 0.001 # converting in Pa*s

## Coolant thermal conductivity - input in K - output in W/m/K
cool_th_cond = 110 - 0.0648*temp + 1.16e-5*temp**2

## Nusselt
Nu_cool = 7 + 0.025*Pe_cool**0.8



#### ************************************ CLADDING ********************************* ####
clad_temp_melting = 1673 # K
clad_temp_max = 650 + 273.15 # K

## Linear thermal expansion: ATTENZIONE usare Kelvin!
clad_eps_th = -3.101e-4 + 1.525e-5*(temp-273.15) + 2.75e-9*(temp-273.15)**2
alfa_clad = 1e-5 # C°^-1 @ 298.15 K

## Cladding density: ATTENZIONE usare Kelvin (see above) - output in kg/m^3
clad_density = 7900*(1+ clad_eps_th )**-3

## Cladding thermal conductivity: Always Kelvin...
clad_thermal_cond = 13.95 + 0.01163*(temp-273.15)

# Clad mech prop
def clad_Young_modulus(temperature):
    #   input: temperature [K]
    #   output: [Pa]

    E=(202.7-0.08167*(temperature-273.15))*1e9
    return E

def clad_Poisson_ratio(temperature):
    #   input: temperature [K]
    #   output: [/]

    nu=0.277+6e-5*(temperature-273.15)
    return nu

#### HELIUM PROPERTIES ####
helium_thermal_cond = 15.8e-4 * temp**0.79


# constants for closed gap
fg_C = 10*(0.0348**-0.5) # m^-0.5
fuel_hardness = 11.2e9 # Pa



#### **************** FUEL PROPERTIES ******************** ####
fuel_temp_max_suggested = 2600 + 273.15 # K

## thermal conductivity: kelvin...
A = 0.01926 + 1.06e-6 * x_om + 2.63e-8 * pu_conc
B = 2.39e-4 + 1.37e-13 * pu_conc
D = 5.27e9
E = 17109.5
k_0 = ( 1/(A + B*temp) + (D/(temp**2))*exp(-E/temp) )*(1-por)**2.5
fuel_thermal_cond = 1.755 + (k_0 - 1.755)*exp( -b_up/128.75 )

## melting temp of fuel pin
fuel_temp_melt = 2964.92 + ( (3147 - 364.85*pu_conc - 1014.15*x_om) - 2964.92 )*exp( -b_up/41.01 )
def fuel_temp_melting(pu=0.29,xom=0,burnup=1e4):
    out = fuel_temp_melt.subs(pu_conc,pu)
    out = out.subs(x_om,xom)
    out = out.subs(b_up,burnup)
    return out

## linear thermal ref fuel
alfa_fuel = 1e-5 # @ 298.15 K

## density
fuel_density = 11.31e3 * 0.945 # kg/m^3

#fuel_E = 250e9 # Pa
fuel_nu = 0.32

def fuel_Young_modulus(temperature,porosity):
    #   input: temperature [K]
    #   porosity [/]
    #   output: [Pa]
    E=((22.43e4-31.19*(temperature-273.15))*(1-2.6*porosity))*1e6
    return E


## restructuring properties
fuel_temp_clmn = 1800 + 273.15 # K
fuel_temp_eqax = 1600 + 273.15 # K

poro_asf = 0.12
poro_clmn = 0.05
poro_void = 1

# (HP SEMPL) considering maximum burnup (i.e. at midplane)
def k_th_fuel(temperature,bup,x=0,pu=0.29,po=0.12):
    """
    BURN UP IN GWd/ton
    """
    output = fuel_thermal_cond.subs(x_om, x)
    output = output.subs(pu_conc, pu)
    output = output.subs(por, po)
    output = output.subs(b_up, bup)
    output = output.subs(temp,temperature)
    return output



def swelling_fuel(z,burnup,size):
    """
    BURN UP IN GWd/ton
    """
    #bup = peak_factor_calc(z)*burnup
    bup = interpolated_peak_factor(z) * burnup

    return (1 + 0.0007/3*bup)*size


def k_th_gas(temperature,x_he=1,x_xe=0,x_kr=0):
    """
    Actually dependance on Burnup (hence time) expressed by mol of gas (x_he, x_xe, x_kr)

    NOTE:
    - BURN UP IN GWd/ton
    """
    x_he_tot = 0.004 + x_he # 4e-4 is an estimation of initial amount of moles
    x_tot = x_he_tot + x_xe + x_kr
    x_he_rel = x_he_tot/x_tot
    x_xe_rel = x_xe/x_tot
    x_kr_rel = x_kr/x_tot

    k_he = 15.8*1e-4 * temperature**0.79
    k_xe = 0.72*1e-4 * temperature**0.79
    k_kr = 1.15*1e-4 * temperature**0.79

    out = k_he**x_he_rel * k_xe**x_xe_rel * k_kr**x_kr_rel
    return out


def swelling_clad(z,burnup,temperature,diam_clad):
    """
    WARNING: output as diameter rather than radius

    NOTE:
    - BURN UP IN GWd/ton
    """
    celsTemp = temperature - 273.15

    #flux = peak_factor_calc(z) * 6.1e15 # n/cm^2/sec
    flux = interpolated_peak_factor(z) * 6.1e15 # n/cm^2/sec

    time = (3600*24)*365 * (burnup/52) # sec - conversion from b-up to seconds (HP constant flux)
    #sw = 1.5e-3 * exp( -2.5 * ( (celsTemp - 450)/100 )**2 ) * ( flux*time/1e22 )**2.75 # %
    sw = 1.3e-5 * exp( -( (celsTemp-490)/100 )**2 ) * ( flux*time/1e22 )**3.9 # USING THE BEST ESTIMATION (SEEN AT LESSON)
    radius = diam_clad/2
    return 2*radius * (1 + 0.01*sw/3)


def yield_strength_cladding(T):
    """
        input: T [°C]
    """
    if T < 600:
        return 555.5e6 - 0.25 * T
    elif T <= 1000:
        return 405.5e6 - 0.775 * (T - 600)
    else:
        return 345.5e6 - 0.25 * T

def UTS_cladding(T):
    """
        input: T [°C]
    """


    if T < 600:
        return 700.5e6 - 0.3215 * T
    elif T <= 1000:
        return 512.5e6 - 0.969 * (T - 600)
    else:
        return 437.5e6 - 0.3125 * T


## FUNCTION USEFUL ABOVE
def peak_factor_calc(z):
    """
    Getting peak factor according to position along the pin
    """
    unit = pin_top_pos / 10  #0.085
    peak_factor = np.array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498])
    value = int(z/unit)

    return peak_factor[value]

def interpolated_peak_factor(z):
    unit = pin_top_pos / 10  #0.085
    peak_factor = np.array([.572, .737, .868, .958, 1, .983, .912, .802, .658, .498, .498]) # ultimo ripetuto poicè si riferisce a z = 0.85
    discrete_power = peak_factor * 1

    discrete_domain = np.arange(pin_bottom_pos + unit / 2, pin_top_pos - unit / 2, unit)
    yy_power = np.zeros_like(discrete_domain)
    for i in range(0, len(discrete_domain)):
        value = int(discrete_domain[i] / unit)
        yy_power[i] = discrete_power[value]

    coeff = np.polyfit(discrete_domain, yy_power, 9)
    poly = np.poly1d(coeff)
    y_new = poly(z)
    return y_new


