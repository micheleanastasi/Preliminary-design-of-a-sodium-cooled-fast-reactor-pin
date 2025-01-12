"""
*** HOT GEOMETRY THERMAL ANALYSIS ***

Script to test and exploit functions and then computing results

"""
import matplotlib.pyplot as plt
import pandas as pd
from functions.thermal_functions import *


burnup = 128/2
res = 11

#### ********************************************* DOMAIN DISCRETIZATION ****************************************** ####


xx = np.linspace(pin_bottom_pos, pin_top_pos, res)
rr = np.linspace(0,fuel_d_outer/2,10)
yy_power_linear = np.zeros_like(xx)
yy_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
yy_hot_temp = np.zeros([len(xx),5]) # the same
yy_gap = np.zeros([len(xx),1])
yy_properties = np.zeros([len(xx),11])
rr_temp_fuel_radial = np.zeros((len(xx),len(rr)))
clad_diam_out = np.zeros_like(xx)
fuel_diam_outer = np.zeros_like(xx)
mean_temp_gap = np.zeros([len(xx),2])
contactPressure = np.zeros_like(xx)


#### **************************************************** CALCS *************************************************** ####

for i in range(0,len(xx)): # Z axis
    yy_power_linear[i] = power_lin_distribution(xx[i])
    yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i], yy_properties[i,:], clad_diam_out[i], fuel_diam_outer[i],contactPressure[i] = hot_geometry_general(
        xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0,burnup)
    mean_temp_gap[i,0] = yy_hot_temp[i,2]
    mean_temp_gap[i,1] = yy_hot_temp[i,3]
   # for j in range(0,len(rr)):  # radius
   #     rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j],xx[i],clad_diam_out[i],fuel_diam_outer[i],clad_thickness_0,burnup)
   #     test += 1

#print(gap_vol_cold())
vol_hot = gap_vol_hot(fuel_diam_outer, clad_diam_out,yy_gap)
#print(vol_hot)
test_1, test_2 = pressure_gap_calc( vol_hot, mean_temp_gap, burnup, plenum_vol=0.85*0.25*pi*clad_d_inner**2,
                                    plenum_clad_d_in=clad_d_inner, temp_plenum=yy_hot_temp[0, 0],print_stuff=True )


## calcolo pressione in funzione di extra volume (in termini di lunghezza)
vol_extra = np.arange(0.5e-6,55e-6,0.5e-6)
res_y = np.zeros_like(vol_extra)
res_x = np.zeros_like(vol_extra)

c = 0
for i in vol_extra:
    res_y[c], res_x[c] = pressure_gap_calc(vol_hot, mean_temp_gap, burnup, plenum_vol=i, plenum_clad_d_in=clad_d_inner,
                                           temp_plenum=yy_hot_temp[0, 0],print_stuff=False)
    c += 1
#print(f"temp di calcolo per extra vol in K: {yy_hot_temp[0,0]}")
plt.figure()
plt.plot(res_x*1000,res_y/1e6)
plt.grid()
plt.xlabel("Lunghezza extra in mm")
plt.ylabel("Pressione totale in cladding")
plt.show()




#### ************************************************* EXCEL PRINT ************************************************ ####
data_hotGeo_tempZ_old = np.array([xx, yy_power_linear,
                                    yy_cold_temp[:,0],yy_cold_temp[:,1],yy_cold_temp[:,2],yy_cold_temp[:,3],yy_cold_temp[:,4]]).T
titles_power_old = ['Position in [m]', 'Linear power [W/m]',
                     'Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]']

data_hotGeo_tempZ_new = np.array([xx, yy_power_linear,
                                   yy_gap[:,0],yy_hot_temp[:,0],yy_hot_temp[:,1],yy_hot_temp[:,2],yy_hot_temp[:,3],yy_hot_temp[:,4],
                                   yy_properties[:,0],yy_properties[:,1],yy_properties[:,2],yy_properties[:,3],yy_properties[:,4],
                                   yy_properties[:,5],yy_properties[:,6],yy_properties[:,7],yy_properties[:,8],yy_properties[:,9],
                                   yy_properties[:,10]]).T
titles_power_new = ['Position in [m]', 'Linear power [W/m]',
                    'Gap in [m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]',
                    'HTC [W/m^2/K]','Re','Pr','Pe','Nu',
                    'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]','spec_heat [J/kg/K]',
                    'th_cond [W/m/K]']

df_power = pd.DataFrame(data_hotGeo_tempZ_old, columns=titles_power_old)
df_power.to_excel(os.path.join("excels","data_hotGeo_alongZ_old.xlsx"),index=False)
df_power = pd.DataFrame(data_hotGeo_tempZ_new, columns=titles_power_new)
df_power.to_excel(os.path.join("excels","data_hotGeo_alongZ_new.xlsx"),index=False)



#### ************************************************* PLOTTING *************************************************** ####


#### ***************** PLOT TEMPERATURES (AXIAL) ******************* ####

## plot coolant, cladding in and out temp for hot and cold geometries
plt.figure()
plt.plot(xx,yy_cold_temp[:,1], label='COLD Cladding external',color='red', linestyle='--')
plt.plot(xx,yy_cold_temp[:,2], label='COLD Cladding internal',color='orange', linestyle='--')
plt.plot(xx,yy_hot_temp[:,0], label='Coolant',color='blue')
plt.plot(xx,yy_hot_temp[:,1], label='HOT Cladding external',color='red')
plt.plot(xx,yy_hot_temp[:,2], label='HOT Cladding internal',color='orange')
plt.plot(xx,np.ones(len(xx))*clad_temp_max, label='Max suggested cladding temp', color='black', linestyle='--')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temp profile of coolant and cladding (inner and outer)")
plt.legend()
plt.grid()
plt.show()

## plot fuel internal and external for cold and hot geometries
plt.figure()
#plt.plot(xx,yy_cold_temp[:,3], label='COLD Fuel external',color='blue', linestyle='--')
#plt.plot(xx,yy_cold_temp[:,4], label='COLD Fuel internal',color='red', linestyle='--')
plt.plot(xx,yy_hot_temp[:,3], label='Fuel external',color='blue')
plt.plot(xx,yy_hot_temp[:,4], label='Fuel internal',color='red')
plt.plot(xx,np.ones(len(xx))*fuel_temp_max_suggested, label='Max suggested fuel temp', color='black', linestyle='--')
plt.plot(xx,np.ones(len(xx))*fuel_temp_melting(burnup=burnup), label='Melting point of fuel (HP CONS))', color='black', linestyle='--')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temp profile of fuel pellet (inner and outer)")
plt.legend()
plt.grid()
plt.show()

## contact pressure
plt.figure(19,figsize=(16, 9))

plt.plot(xx,contactPressure*1e-6,label='Contact p @ 52 GWd/ton')

plt.xlabel("Position in [m]")
plt.ylabel("Pressure in [MPa]")
plt.title("Contact pressure w.r.t. burnup")
plt.legend(loc='best')
plt.grid()
plt.show()


## delta tra clad in e out (hot)
plt.figure()
delta_temp_cl = yy_hot_temp[:,2] - yy_hot_temp[:,1]
plt.plot(xx,delta_temp_cl)
plt.xlabel("Position in [m]")
plt.ylabel("Delta temperature in [K]")
plt.title("Axial delta temp btw clad in e clad out @ HOT GEO")
plt.grid()
plt.show()

#### ***************** PLOT TEMPERATURES (RADIAL) ******************* ####
plt.figure()
plt.plot(rr,rr_temp_fuel_radial[int( len(xx)/2 ),:], label="Temp profile")
plt.plot(rr,np.ones(len(rr))*fuel_temp_max_suggested, label="Max suggested fuel temp", color='black', linestyle='--')
plt.plot(rr,np.ones(len(rr))*fuel_temp_melting(burnup=burnup), label="Max suggested fuel temp", color='black', linestyle='--')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Middle position of pin, fuel pellet temperature profile")
plt.grid()
plt.show()

## plot 3d temp radial profile of fuel
x,y = np.meshgrid(rr,xx)
fig_1 = plt.figure()
ax = fig_1.add_subplot(111, projection='3d')
ax.plot_surface(x,y,rr_temp_fuel_radial,cmap='viridis')
ax.set_ylabel('Position along the pin [m]')
ax.set_xlabel('Radius [m]')
ax.set_zlabel('Temperature [K]')
plt.show()


#### ***************** OTHER PROPERTIES (AXIAL) ******************* ####
## GAP
before_gap = ( clad_d_outer - 2*clad_thickness_0 - fuel_d_outer )/2
plt.figure()
plt.plot(xx,yy_gap, label='Gap @ hot geometry')
plt.plot(xx,np.ones_like(xx)*before_gap, label='Gap @ cold geometry', linestyle='--', color='black')
plt.plot(xx,np.zeros_like(xx), label='Gap @ cold geometry', linestyle='--', color='black')
plt.xlabel("Position [m]")
plt.ylabel("Gap [m]")
plt.title("GAP along the pin")
plt.grid()
plt.show()

## velocity of coolant along pin
max_vel = 8 # m/s
plt.figure()
plt.plot(xx,yy_properties[:,5], label='Velocity')
plt.plot(xx,np.ones_like(xx)*max_vel, label='Maximum velocity allowed', linestyle='--', color='black')
plt.xlabel("Position [m]")
plt.ylabel("Velocity [m/s]")
plt.title("Average velocity (over the z section), HOT GEOMETRY")
plt.legend()
plt.grid()
plt.show()