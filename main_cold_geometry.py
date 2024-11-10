
import pandas as pd
from thermal_functions import *
import matplotlib.pyplot as plt

#### ******************** DATA ARRAY INITIALIZATION ************************* ####
xx = np.linspace(pin_bottom_pos,pin_top_pos,100)
rr = np.linspace(0,fuel_d_outer/2,1)

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
test = 0

for i in range(0,len(xx)):
    yy_power_linear[i] = power_lin_distribution(xx[i])
    yy_temp_coolant[i] = temp_coolant(xx[i])
    yy_htc_local[i],yy_adim_number_cool[:,i],yy_cool_local_prop[:,i] = heat_transfer_coeff_local(yy_temp_coolant[i],clad_d_outer)
    yy_temp_clad_out[i] = temp_cladding_outer(xx[i],clad_d_outer)
    yy_temp_clad_in[i] = temp_cladding_inner(xx[i],clad_d_outer,clad_thickness_0)
    yy_temp_fuel_out[i],_ = temp_fuel_outer(xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)
    yy_temp_fuel_in[i] = temp_fuel_inner(xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)

    for j in range(0,len(rr)):
        rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j],xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0)
        test += 1
        print(f"Radial temp pellet calc complete at {np.round(100*test/(len(xx)*len(rr)),3)} %")



# *********************** OUTPUTS ************************* #

data_coldGeo_tempZ = np.array([ xx,yy_power_linear,yy_temp_coolant, yy_temp_clad_out, yy_temp_clad_in, yy_temp_fuel_out, yy_temp_fuel_in,yy_htc_local,yy_adim_number_cool[0,:],yy_adim_number_cool[1,:],yy_adim_number_cool[2,:],yy_adim_number_cool[3,:],
                                                yy_cool_local_prop[0,:],yy_cool_local_prop[1,:],yy_cool_local_prop[2,:],yy_cool_local_prop[3,:],yy_cool_local_prop[4,:]]).T
titles_power = ['Position in [m]','Linear power [W/m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]','HTC [W/m^2/K]','Re','Pr','Pe','Nu',
                'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]','spec_heat [J/kg/K],th_cond [W/m/K]']

#data_coldGeo_tempR = np.array([ rr,yy_power_linear,rr_temp_fuel_radial ])
#df_radial = pd.DataFrame(data_coldGeo_tempR)
#df_radial.to_excel('data_coldGeo_tempR.xlsx',index=False)

df_power = pd.DataFrame(data_coldGeo_tempZ, columns=titles_power)
df_power.to_excel("data_coldGeo_stepPower_power_tempZ.xlsx",index=False)


# ************************ PLOTS ************************** #

#### PLOT POWER DISTRIBUTION ####

plt.figure()
plt.plot(xx,yy_power_linear, label='Power distribution')
plt.xlabel("Position in [m]")
plt.ylabel("Linear power density in [W/m]")
plt.legend()
plt.grid()
plt.show()


#### *********************** PLOT TEMPERATURES (AXIAL) ************************** ####


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


#### *********************** PLOT TEMPERATURES (RADIAL) ************************** ####
# NB cold geo, no redistr, no restructuring, no burn up...
plt.figure()
plt.plot(rr,rr_temp_fuel_radial[int( len(xx)/2 ),:])
plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Cold geometry, middle position (pin) fuel pellet temperature")
plt.grid()
plt.show()

x,y = np.meshgrid(xx,rr)
fig_1 = plt.figure()
ax = fig_1.add_subplot(111, projection='3d')
ax.plot_surface(x,y,rr_temp_fuel_radial,cmap='viridis')
ax.set_ylabel('Position along the pin [m]')
ax.set_xlabel('Radius [m]')
ax.set_zlabel('Temperature [K]')
plt.show()