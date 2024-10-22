
import sympy as sy
from sympy import exp

#### COOLANT PROPERTIES FUNCTIONS AND ASSOCIATED CORRELATIONS ####
cool_T_melting_atm = 98 + 273.15 # K @ p = 0.1 Mpa
cool_T_boiling_atm = 882 + 273.15 # K @ p = 0.1 Mpa

## Symbols used - conversioni temperatura
temp = sy.Symbol('T[K]')
temp_f = sy.Symbol('T[°F]')
Pe = sy.Symbol('Pe')

## Coolant specific heat - input in K - output in J/kg/K
cool_spec_heat = 1608 - 0.7481*temp + 3.929e-4*temp ** 2

## Coolant density - input in °F - output in kg/m^3
cool_density__fahr = 954.1579 + temp_f*( temp_f*(temp_f*0.9667e-9 - 0.46e-5) - 0.1273534 )

## Coolant dyn viscosity - input in K - output in Pa*s
cool_dyn_viscosity = exp( 813.9/temp - 2.530 ) * 0.001 # converto in Pa*s

## Coolant thermal conductivity - input in K - output in W/m/K
cool_th_cond = 110 - 0.0648*temp + 1.16e-5*temp**2

## Nusselt - VALIDITà????
Nu = 7 + 0.025*Pe**0.8



#### CLADDING ####
