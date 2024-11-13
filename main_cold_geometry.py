"""
*** COLD GEOMETRY ANALYSIS ***

NOTE: this approach is the most conservative ever, since we're neglecting a lot, if not all, of beneficial phenomena...
"""
import pandas as pd
from thermal_functions import *
import matplotlib.pyplot as plt

#### ******************************************* DOMAIN DISCRETIZATION ******************************************** ####
xx = np.linspace(pin_bottom_pos,pin_top_pos,20)
rr = np.linspace(0,fuel_d_outer/2,10)



#### ******************************************* DATA ARRAY INITIALIZATION **************************************** ####
yy_power_linear = np.zeros_like(xx)
yy_htc_local = np.zeros_like(xx)
yy_adim_number_cool = np.zeros([4,len(xx)])
yy_cool_local_prop = np.zeros([6,len(xx)])
yy_temp_coolant = np.zeros_like(xx)
yy_temp_clad_out = np.zeros_like(xx)
yy_temp_clad_in = np.zeros_like(xx)
yy_temp_fuel_out = np.zeros_like(xx)
yy_temp_fuel_in = np.zeros_like(xx)
rr_temp_fuel_radial = np.zeros((len(xx),len(rr)))
test = 0     # to keep tracking...



#### **************************************************** CALCS *************************************************** ####
for i in range(0,len(xx)):
    yy_power_linear[i] = power_lin_distribution(xx[i])
    yy_temp_coolant[i] = temp_coolant(xx[i])
    yy_htc_local[i],yy_adim_number_cool[:,i],yy_cool_local_prop[:,i] = heat_transfer_coefficient(yy_temp_coolant[i],clad_d_outer)
    yy_temp_clad_out[i] = temp_cladding_outer(xx[i],clad_d_outer)
    yy_temp_clad_in[i] = temp_cladding_inner(xx[i],clad_d_outer,clad_thickness_0)
    yy_temp_fuel_out[i],_ = temp_fuel_outer(xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)
    yy_temp_fuel_in[i] = temp_fuel_inner(xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)

    for j in range(0,len(rr)):
        rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j],xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0)
        test += 1    # to keep tracking...
        print(f"Radial temp pellet calc complete at {np.round(100*test/(len(xx)*len(rr)),3)} %")



#### ******************************************** EXCEL OUTPUTS *************************************************** ####

data_coldGeo_tempZ = np.array([ xx, yy_power_linear, yy_temp_coolant, yy_temp_clad_out,
                                yy_temp_clad_in, yy_temp_fuel_out, yy_temp_fuel_in, yy_htc_local,
                                yy_adim_number_cool[0,:],yy_adim_number_cool[1,:],yy_adim_number_cool[2,:],yy_adim_number_cool[3,:],
                                yy_cool_local_prop[0,:],yy_cool_local_prop[1,:],yy_cool_local_prop[2,:],yy_cool_local_prop[3,:],
                                yy_cool_local_prop[4,:],yy_cool_local_prop[5,:]]).T
titles_power = ['Position in [m]','Linear power [W/m]','Temp coolant [K]','Temp cladding outer [K]',
                'Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]','HTC [W/m^2/K]',
                'Re','Pr','Pe','Nu',
                'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]',
                'spec_heat [J/kg/K]','th_cond [W/m/K]']

df_power = pd.DataFrame(data_coldGeo_tempZ, columns=titles_power)
df_power.to_excel("data_coldGeo_alongZ.xlsx",index=False)



#### ************************************************** PLOTS ***************************************************** ####

#### PLOT POWER DISTRIBUTION ####

plt.figure(1)
plt.plot(xx,yy_power_linear, label='Power distribution')
plt.xlabel("Position in [m]")
plt.ylabel("Linear power density in [W/m]")
plt.legend()
plt.grid()
plt.show()


#### *********************** PLOT TEMPERATURES (AXIAL) ************************** ####


plt.figure(2)
plt.plot(xx,yy_temp_coolant, label='Coolant')
plt.plot(xx,yy_temp_clad_out, label='Cladding external')
plt.plot(xx,yy_temp_clad_in, label='Cladding internal')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.legend()
plt.grid()
plt.show()


plt.figure(3)
plt.plot(xx,yy_temp_fuel_out, label='Fuel external')
plt.plot(xx,yy_temp_fuel_in, label='Fuel internal')
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.legend()
plt.grid()
plt.show()


#### *********************** PLOT TEMPERATURES (RADIAL) ************************** ####
# NB cold geo, no redistr, no restructuring, no burn up...
plt.figure(4)
plt.plot(rr,rr_temp_fuel_radial[int( len(xx)/2 ),:])
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Cold geometry, middle position (pin) fuel pellet temperature")
plt.grid()
plt.show()

# 3d plot
x,y = np.meshgrid(rr,xx)
fig_1 = plt.figure(5)
ax = fig_1.add_subplot(111, projection='3d')
ax.plot_surface(x,y,rr_temp_fuel_radial,cmap='viridis')
ax.set_ylabel('Position along the pin [m]')
ax.set_xlabel('Radius [m]')
ax.set_zlabel('Temperature [K]')
plt.show()