"""
Material properties deriving from homework.pdf data
"""

import sympy as sy
from sympy import exp

## Symbols used from sympy - conversioni temperatura
temp = sy.Symbol('T[K]')
temp_f = sy.Symbol('T[°F]')
Pe_cool = sy.Symbol('Pe')
x_om = sy.Symbol('x[O/M]')
pu_conc = sy.Symbol('[Pu]')
por = sy.Symbol('P')



#### ****************** COOLANT PROPERTIES ********************* ####

cool_T_melting_atm = 98 + 273.15 # K @ p = 0.1 Mpa
cool_T_boiling_atm = 882 + 273.15 # K @ p = 0.1 Mpa

## Coolant specific heat - input in K - output in J/kg/K
cool_spec_heat = 1608 - 0.7481*temp + 3.929e-4*temp ** 2

## Coolant density - input in °F!!!!!! - output in kg/m^3
cool_density__fahr = 954.1579 + temp_f*( temp_f*(temp_f*0.9667e-9 - 0.46e-5) - 0.1273534 )

## Coolant dyn viscosity - input in K - output in Pa*s
cool_dyn_viscosity = exp( 813.9/temp - 2.530 ) * 0.001 # converto in Pa*s

## Coolant thermal conductivity - input in K - output in W/m/K
cool_th_cond = 110 - 0.0648*temp + 1.16e-5*temp**2

## Nusselt - VERIFICA VALIDITà????
Nu_cool = 7 + 0.025*Pe_cool**0.8



#### ********** CLADDING *********** ####
clad_temp_melting = 1673 # K
clad_temp_max = 650 + 273.15 # K

## Linear thermal expansion: ATTENZIONE usare Kelvin!
clad_eps_th = -3.101e-4 + 1.525e-5*(temp-273.15) + 2.75e-9*(temp-273.15)**2
alfa_clad = 1e-5 # C°^-1 @ 298.15 K

## Cladding density: ATTENZIONE usare Kelvin (see above) - output in kg/m^3
clad_density = 7900*(1+ clad_eps_th )**-3

## Cladding thermal conductivity: Always Kelvin...
clad_thermal_cond = 13.95 + 0.01163*(temp-273.15)

#### HELIUM PROPERTIES ####
helium_thermal_cond = 15.8e-4 * temp**0.79



#### **************** FUEL PROPERTIES ******************** ####
fuel_temp_max_suggested = 2600 + 273.15 # K

## thermal conductivity: kelvin...
# per adesso usare x = 2, Pu = 20%, por = 12%
A = 0.01926 + 1.06e-6 * x_om + 2.63e-8 * pu_conc
B = 2.39e-4 + 1.37e-13 * pu_conc
D = 5.27e9
E = 17109.5
k_0 = ( 1/(A + B*temp) + (D/(temp**2))*exp(-E/temp) )*(1-por)**2.5

# da aggiungere burn up dopo! ( al posto di 0 --> hp conservativa???? )
fuel_thermal_cond = 1.755 + (k_0 - 1.755)*0

## melting temp
# da aggiungere burn up dopo! ( al posto di 0 --> hp conservativa???? )
fuel_temp_melting = 2964.92 + ( (3147 - 364.85*pu_conc - 1014.15*x_om) - 2964.92 )*0

## linear thermal ref fuel
alfa_fuel = 1e-5 # @ 298.15 K

## restructuring temps
fuel_temp_clmn = 1800 + 273.15 # K
fuel_temp_eqax = 1600 + 273.15 # K

