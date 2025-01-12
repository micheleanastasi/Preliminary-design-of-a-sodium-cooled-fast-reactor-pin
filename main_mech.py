"""
GENERAL MAIN SCRIPT USED TO CALCULATE EVERYTHING CONCERNING MECHANICAL ASPECTS!
"""
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import pandas as pd
from functions.general_functions import *
from functions.general_properties import *
from functions.mechanical_functions import *
import os

yy_power_linear = np.load(os.path.join("main_numpy_saves", "yy_power_linear.npy"))  # linear power
yy_hot_temp = np.load(os.path.join("main_numpy_saves", "yy_hot_temp.npy"))  #   axial temperature
yy_cold_temp = np.load(os.path.join("main_numpy_saves", "yy_cold_temp.npy"))    #   coolant temperature
yy_gap = np.load(os.path.join("main_numpy_saves", "yy_gap.npy"))    #   gap size
yy_properties = np.load(os.path.join("main_numpy_saves", "yy_properties.npy"))  #   'other' in hot_geometry_general: heat transfer coeff, adimensional numbers, coolant properties
clad_diam_out = np.load(os.path.join("main_numpy_saves", "clad_diam_out.npy"))
fuel_diam_outer = np.load(os.path.join("main_numpy_saves", "fuel_diam_outer.npy"))
vol_hot = np.load(os.path.join("main_numpy_saves", "vol_hot.npy"))
pressure = np.load(os.path.join("main_numpy_saves", "pressure.npy"))
extra_pin = np.load(os.path.join("main_numpy_saves", "extra_pin.npy"))
yy_r_clmn = np.load(os.path.join("main_numpy_saves", "yy_r_clmn.npy"))
yy_r_void = np.load(os.path.join("main_numpy_saves", "yy_r_void.npy"))
sw_clad = np.load(os.path.join("main_numpy_saves", "sw_clad.npy"))
contactPress = np.load(os.path.join("main_numpy_saves", "contactPress.npy"))


burnup=64   #[Gwd/t]


## VARIABLE INITIALIZATION

yy_temp_clad_in=yy_hot_temp[:,2,2]
stress_yield=np.zeros(yy_temp_clad_in.shape)
stress_rupture=np.zeros(yy_temp_clad_in.shape)
stress_yield_embrittled=np.zeros(yy_temp_clad_in.shape)
stress_rupture_embrittled=np.zeros(yy_temp_clad_in.shape)

helium_content_moles = helium_content(burnup)   #He content in cladding


for i in range(len(yy_temp_clad_in)):
    stress_yield[i]=yield_strength_cladding(yy_temp_clad_in[i]-273.15)
    stress_yield_embrittled[i] = yield_stress_embrittlement(stress_yield[i], helium_content_moles)
    stress_rupture[i]=UTS_cladding(yy_temp_clad_in[i]-273.15)
    stress_rupture_embrittled[i] = yield_stress_embrittlement(stress_rupture[i], helium_content_moles)

print(f'Yield stress @ Bu=64 [GWd/t] as a function of z : {stress_yield_embrittled*1e-6} [Mpa]')
print(f'Rupture stress at Bu = 64 [GWd/t] : {max(stress_rupture_embrittled)*1e-6} [MPa]')

gas_pressure=pressure[2]
total_pressure=contactPress[:,2]+gas_pressure
avg_pressure=np.average(total_pressure)


print('******************')
print('Contact pressure @ Bu=64 [GWd/t] \\')
print(f'Gas pressure = {gas_pressure*1e-6} [MPa]')
print(f'Total (avg) pressure = {avg_pressure*1e-6} [MPa]')
print(f'Max pressure = {max(total_pressure)*1e-6} [MPa]')



#### ***************** cladding thickness - Mariotte solution - Tresca criterion

stress_r=-avg_pressure/2
stress_theta=avg_pressure*clad_diam_out[:,2]/2/clad_thickness_0
stress_z=avg_pressure*clad_diam_out[:,2]/2/2/clad_thickness_0

# critical_stress=avg_pressure*(clad_diam_out[:,2]/2/clad_thickness_0+0.5)
critical_stress=np.abs(stress_theta-stress_r)
print(f'Critical stress  = {max(critical_stress)*1e-6} [MPa]')

max_total_pressure=stress_yield_embrittled/(clad_diam_out[:,2]/2/clad_thickness_0+0.5)
critical_stress*1e-6


##  Thermal stresses

clad_E=clad_Young_modulus((yy_hot_temp[:,2,2]+yy_hot_temp[:,1,2])/2)
clad_nu=clad_Poisson_ratio((yy_hot_temp[:,2,2]+yy_hot_temp[:,1,2])/2)
clad_r_out=clad_d_outer/2
clad_diam_in=clad_d_outer-2*clad_thickness_0
clad_r_in=clad_diam_out/2
clad_r_avg=(clad_r_in+clad_r_out)/2

th_stress=alfa_clad*clad_E/(1-clad_nu)*(max(yy_hot_temp[:,2,2]-yy_hot_temp[:,1,2])/2)

print(f'Thermal stress as a function of z : {(th_stress)*1e-6} [MPa]')


## Time to rupture

crit_index=np.where(critical_stress==max(critical_stress))[0][0]
crit_T=yy_temp_clad_in[crit_index]
LMP=(2060-max(critical_stress+max(th_stress))*1e-6)/0.095

time_to_rupture=10**(LMP/crit_T-17.125)

print(f'Time to rupture = {time_to_rupture/24} days')