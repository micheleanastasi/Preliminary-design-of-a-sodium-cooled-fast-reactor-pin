import numpy as np
import pandas as pd
from thermal_functions import *

xx = np.linspace(pin_bottom_pos,pin_top_pos,100)

#### DATA GUESS ####
clad_thickness_0 = 0.5e-3 # m
clad_fuel_gap = 80e-6 # m


#### DATA ARRAY ####
yy_power_linear = np.zeros_like(xx)

yy_temp_coolant = np.zeros_like(xx)
yy_temp_clad_out = np.zeros_like(xx)
yy_temp_clad_in = np.zeros_like(xx)
yy_temp_fuel_out = np.zeros_like(xx)
yy_temp_fuel_in = np.zeros_like(xx)


for i in range(0,len(xx)):
    yy_power_linear[i] = power_lin_distribution(xx[i])

    yy_temp_coolant[i] = temp_coolant(xx[i])
    yy_temp_clad_out[i] = temp_cladding_outer(xx[i],clad_d_outer)
    yy_temp_clad_in[i] = temp_cladding_inner(xx[i],clad_d_outer,clad_thickness_0)
    yy_temp_fuel_out[i] = temp_fuel_outer(xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0,clad_fuel_gap)
    yy_temp_fuel_in[i] = temp_fuel_inner(xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0,clad_fuel_gap)

# ********************************************************* #
# *********************** OUTPUTS ************************* #

data_coldGeo_stepPower_power_tempZ = np.array([ xx,yy_power_linear,yy_temp_coolant, yy_temp_clad_out, yy_temp_clad_in, yy_temp_fuel_out, yy_temp_fuel_in ]).T
titles_power = ['Position in [m]','Linear power [W/m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]']

df_power = pd.DataFrame(data_coldGeo_stepPower_power_tempZ, columns=titles_power)
df_power.to_excel("data_coldGeo_stepPower_power_tempZ.xlsx",index=False)


# ********************************************************* #
# ************************ PLOTS ************************** #

#### PLOT POWER DISTRIBUTION ####

plt.figure()
plt.plot(xx,yy_power_linear, label='Power distribution')
plt.xlabel("Position in [m]")
plt.ylabel("Linear power density in [W/m]")
plt.legend()
plt.grid()
plt.show()


#### PLOT TEMPERATURES (AXIAL) ####


plt.figure()
plt.plot(xx,yy_temp_coolant, label='Coolant')
plt.plot(xx,yy_temp_clad_out, label='Cladding external')
plt.plot(xx,yy_temp_clad_in, label='Cladding internal')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.legend()
plt.grid()
plt.show()


plt.figure()
plt.plot(xx,yy_temp_fuel_out, label='Fuel external')
plt.plot(xx,yy_temp_fuel_in, label='Fuel internal')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.legend()
plt.grid()
plt.show()