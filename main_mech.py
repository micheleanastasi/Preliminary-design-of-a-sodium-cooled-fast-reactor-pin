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


burnup=52   #[Gwd/t]


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

print(f'Yield stress @ Bu=52 [GWd/t] : {stress_yield_embrittled*1e-6} [Mpa]')


print(f'Max contact pressure at Bu = 52 [GWd/t] : {max(contactPress[:,2])*1e-6} [MPa]')
gas_pressure=pressure[2]
total_pressure=contactPress[:,2]+gas_pressure
print(f'Total pressure at Bu= 52 [GWd/t] : {max(total_pressure)*1e-6} [MPa]')



#### ***************** cladding thickness - Mariotte solution - Tresca criterion

critical_stress=total_pressure*(clad_diam_out[:,2]/2/clad_thickness_0+0.5)
print(f'Critical stress at Bu = 52 [GWd/t] : {max(critical_stress)*1e-6} [MPa]')

max_total_pressure=stress_yield_embrittled/(clad_diam_out[:,2]/2/clad_thickness_0+0.5)
critical_stress*1e-6
# 120*(clad_diam_out[:,2]/2/clad_thickness_0+0.5)

stress_rupture_embrittled*1e-6

## Time to rupture

crit_index=np.where(critical_stress==max(critical_stress))[0][0]
crit_T=yy_temp_clad_in[crit_index]
LMP=(2060-max(critical_stress)*1e-6)/0.095

time_to_rupture=10**(LMP/crit_T-17.125)