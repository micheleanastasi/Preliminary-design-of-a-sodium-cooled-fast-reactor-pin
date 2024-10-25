import numpy as np
import pandas as pd
from PIL.ImageColor import colormap

from thermal_functions import *

xx = np.linspace(pin_bottom_pos,pin_top_pos,100)



def hot_geometry_iteration(z,clad_d_out_0,fuel_d_out_0,clad_thick_0):
    tol = 10e-3
    temp_array = np.zeros(5)
    delta_gap = (clad_d_out_0 - 2 * clad_thick_0 - fuel_d_out_0) / 2

    temp_array[0] = temp_coolant(z)
    temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
    temp_array[2] = temp_cladding_inner(z, clad_d_out_0, delta_gap)
    temp_array[3] = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
    temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
    old = temp_array.copy()

    while True:
        prec_clad_d_out = clad_d_out_0
        prec_fuel_d_out = fuel_d_out_0
        prec_temp_array = temp_array.copy()
        prec_delta_gap = delta_gap

        clad_d_out_0 = diameter_th_exp_cladding(clad_d_outer, prec_temp_array[1]) # temp clad outer - DIAMETRO INIZIALE!!!
        fuel_d_out_0 = diameter_th_exp_fuel(fuel_d_outer, prec_temp_array[3]) # temp fuel outer - DIAMETRO INIZIALE !!!

        delta_gap = (clad_d_out_0 - 2 * clad_thick_0 - fuel_d_out_0) / 2

        temp_array[0] = temp_coolant(z)
        temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
        temp_array[2] = temp_cladding_inner(z, clad_d_out_0, delta_gap)
        temp_array[3] = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
        temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)

        if np.abs(prec_temp_array[4] - temp_array[4] ) < tol :
            break
    print(f"Completed at {np.round(100*z/0.85,2)}% (Position: {np.round(z,2)} m) - Temp fuel inner: HOT:{np.round(temp_array[4],2)}, COLD:{np.round(old[4],2)} - Gap:{np.round(delta_gap*1000,6)} mm")
    return old,temp_array,delta_gap

#### ******* OUTPUT ******* ####

yy_cold_temp = np.zeros([5,len(xx)])
yy_hot_temp = np.zeros([5,len(xx)])
yy_gap = np.zeros([1,len(xx)])

for i in range(0,len(xx)):
    yy_cold_temp[:,i], yy_hot_temp[:,i],yy_gap = hot_geometry_iteration(xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0)

#data_hotGeo_stepPower_tempZ_cold = np.array([xx,yy_cold_temp[0],yy_cold_temp[1],yy_cold_temp[2],yy_cold_temp[3],yy_cold_temp[4]]).T
#data_hotGeo_stepPower_tempZ_hot = np.array([xx,yy_gap,yy_hot_temp[0],yy_hot_temp[1],yy_hot_temp[2],yy_hot_temp[3],yy_hot_temp[4]]).T

#titles_power_cold = ['Position in [m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]']
#titles_power_hot = ['Position in [m]','Gap in [m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]']

#df_power = pd.DataFrame(data_hotGeo_stepPower_tempZ_cold, columns=titles_power_cold)
#df_power.to_excel("data_hotGeo_stepPower_tempZ_cold.xlsx",index=False)

#df_power = pd.DataFrame(data_hotGeo_stepPower_tempZ_hot, columns=titles_power_hot)
#df_power.to_excel("data_hotGeo_stepPower_tempZ_hot.xlsx",index=False)


#### PLOT TEMPERATURES (AXIAL) ####


plt.figure()
plt.plot(xx,yy_cold_temp[3], label='COLD Fuel external',color='blue', linestyle='--')
plt.plot(xx,yy_cold_temp[4], label='COLD Fuel internal',color='red', linestyle='--')
plt.plot(xx,yy_hot_temp[3], label='HOT Fuel external',color='blue')
plt.plot(xx,yy_hot_temp[4], label='HOT Fuel internal',color='red')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Cold and hot geometry")
plt.legend()
plt.grid()
plt.show()