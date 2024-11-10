"""
*** HOT GEOMETRY THERMAL ANALYSIS ***

NOTE: still missing phenomena such as restructuring, burn-up and redistribution (the latter has minor impact though)
    Nonetheless, up to this point everything was taken into account in a conservative approach, although some aspect should
    be evaluated accordingly (i.e. how much is Pu redistr. negligible?)
    But why conservative? Burn up -> inf, then lowering the values of T_fuel_melting, k_thermal, etc at minimum... moreover not still
    considering restructuring!
    THIS IS THE FIRST STEP

SO:
- start: only hot geometry, no burn up (we are HERE but cons.hp: burn up --> inf)
- then: hot geo + restructuring + Pu redistribution
- end: hot geo + "" + "" + burn up
properties changing: k_fuel, Tm_fuel, ...
"""
import matplotlib.pyplot as plt
import pandas as pd
from thermal_functions import *


#### ***************** DOMAIN DISCRETIZATION ********************** ####
resol = 20
xx = np.linspace(pin_bottom_pos,pin_top_pos,resol)
rr = np.linspace(0,fuel_d_outer/2,resol)



#### *************************** HOT GEOMETRY CALCULATIONS ************************** ####
def hot_geometry_iteration(z,clad_d_out_0,fuel_d_out_0,clad_thick_0):
    tol = 10e-3

    # initialize variables
    temp_array = np.zeros(5)
    delta_gap = (clad_d_out_0 - 2*clad_thick_0 - fuel_d_out_0)/2
    yy_htc_loc = np.zeros([1])
    yy_adim_num_cool = np.zeros([4])
    yy_cool_loc_prop = np.zeros([6])

    temp_array[0] = temp_coolant(z)
    temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
    yy_htc_loc, yy_adim_num_cool, yy_cool_loc_prop = heat_transfer_coeff_local(temp_array[0],clad_d_out_0) # new
    temp_array[2] = temp_cladding_inner(z, clad_d_out_0, delta_gap)
    temp_array[3],_ = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
    temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
    old = temp_array.copy()

    while True:
        prec_clad_d_out = clad_d_out_0
        prec_fuel_d_out = fuel_d_out_0
        prec_temp_array = temp_array.copy()
        prec_delta_gap = delta_gap

        clad_d_out_0 = diameter_th_exp_cladding(clad_d_outer, prec_temp_array[1]) # temp clad outer - DIAMETRO INIZIALE!!!
        fuel_d_out_0 = diameter_th_exp_fuel(fuel_d_outer, prec_temp_array[3]) # temp fuel outer - DIAMETRO INIZIALE !!!

        temp_array[0] = temp_coolant(z)
        temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
        yy_htc_loc, yy_adim_num_cool, yy_cool_loc_prop = heat_transfer_coeff_local(temp_array[0],clad_d_out_0)  # new
        temp_array[2] = temp_cladding_inner(z, clad_d_out_0, delta_gap)
        temp_array[3],delta_gap = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)
        temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thickness_0)

        if np.abs(prec_temp_array[4] - temp_array[4] ) < tol : # va bene così (?)
            break

    print(f"Hot geo completed at {np.round(100*z/0.85,2)}% (Position: {np.round(z,2)} m) - Temp fuel inner: HOT:{np.round(temp_array[4],2)}, COLD:{np.round(old[4],2)} - Gap:{delta_gap*1000} mm")
    other = np.array( list([yy_htc_loc]) + list(yy_adim_num_cool) + list(yy_cool_loc_prop) )
    return old,temp_array,delta_gap,other, clad_d_out_0, fuel_d_out_0



#### *********************** CALCS *********************** ####

yy_power_linear = np.zeros_like(xx)
yy_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
yy_hot_temp = np.zeros([len(xx),5]) # the same
yy_gap = np.zeros([len(xx),1])
yy_properties = np.zeros([len(xx),11])
rr_temp_fuel_radial = np.zeros((len(xx),len(rr)))

test = 0
for i in range(0,len(xx)):
    yy_power_linear[i] = power_lin_distribution(xx[i])
    yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i],yy_properties[i,:],clad_diam_out,fuel_diam_outer = hot_geometry_iteration(xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0)
    for j in range(0,len(rr)):
        rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j],xx[i],clad_diam_out,fuel_diam_outer,clad_thickness_0)
        test += 1
        #print(f"Calcs completed at {np.round(100*test/(len(xx)*len(rr)),3)} %")



#### ********************** EXCEL PRINT ********************* ####

data_hotGeo_tempZ_cold = np.array([ xx,yy_power_linear,yy_cold_temp[:,0],yy_cold_temp[:,1],yy_cold_temp[:,2],yy_cold_temp[:,3],yy_cold_temp[:,4],yy_properties[:,0],
                                              yy_properties[:,1],yy_properties[:,2],yy_properties[:,3],yy_properties[:,4],yy_properties[:,5],yy_properties[:,6],yy_properties[:,7],yy_properties[:,8],yy_properties[:,9],yy_properties[:,10] ]).T
data_hotGeo_tempZ_hot = np.array([ xx,yy_power_linear,yy_gap[:,0],yy_hot_temp[:,0],yy_hot_temp[:,1],yy_hot_temp[:,2],yy_hot_temp[:,3],yy_hot_temp[:,4],yy_properties[:,0],
                                              yy_properties[:,1],yy_properties[:,2],yy_properties[:,3],yy_properties[:,4],yy_properties[:,5],yy_properties[:,6],yy_properties[:,7],yy_properties[:,8],yy_properties[:,9],yy_properties[:,10] ]).T
titles_power_cold = ['Position in [m]','Linear power [W/m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]','HTC [W/m^2/K]','Re','Pr','Pe','Nu',
                'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]','spec_heat [J/kg/K]','th_cond [W/m/K]']
titles_power_hot = ['Position in [m]','Linear power [W/m]','Gap in [m]','Temp coolant [K]','Temp cladding outer [K]','Temp cladding inner [K]','Temp fuel outer [K]','Temp fuel inner [K]','HTC [W/m^2/K]','Re','Pr','Pe','Nu',
                'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]','spec_heat [J/kg/K]','th_cond [W/m/K]']

df_power = pd.DataFrame(data_hotGeo_tempZ_cold, columns=titles_power_cold)
df_power.to_excel("data_hotGeo_interPower_tempZ_cold.xlsx",index=False)
df_power = pd.DataFrame(data_hotGeo_tempZ_hot, columns=titles_power_hot)
df_power.to_excel("data_hotGeo_interPower_tempZ_hot.xlsx",index=False)


def plotting():
    #### ***************** PLOT TEMPERATURES (AXIAL) ******************* ####
    ## NB no redistr, no restructuring, no burn up...

    # plot coolant, cladding in and out temp for hot and cold geometries
    plt.figure()
    #plt.plot(xx,yy_cold_temp[:,0], label='COLD Coolant',color='blue', linestyle='--')
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

    # plot fuel internal and external for cold and hot geometries
    plt.figure()
    plt.plot(xx,yy_cold_temp[:,3], label='COLD Fuel external',color='blue', linestyle='--')
    plt.plot(xx,yy_cold_temp[:,4], label='COLD Fuel internal',color='red', linestyle='--')
    plt.plot(xx,yy_hot_temp[:,3], label='HOT Fuel external',color='blue')
    plt.plot(xx,yy_hot_temp[:,4], label='HOT Fuel internal',color='red')
    plt.plot(xx,np.ones(len(xx))*fuel_temp_max_suggested, label='Max suggested fuel temp', color='black', linestyle='--')
    plt.plot(xx,np.ones(len(xx))*2964.92, label='Melting point of fuel (HP CONS))', color='black', linestyle='--')
    plt.xlabel("Position in [m]")
    plt.ylabel("Temperature in [K]")
    plt.title("Axial temp profile of fuel pellet (inner and outer)")
    plt.legend()
    plt.grid()
    plt.show()


    # plot difference btw fuel hot and cold...
    plt.figure()
    plt.plot(xx,yy_cold_temp[:,3] - yy_hot_temp[:,3], label='Fuel external',color='blue')
    plt.plot(xx,yy_cold_temp[:,4] - yy_hot_temp[:,4], label='Fuel internal',color='red')
    plt.xlabel("Position in [m]")
    plt.ylabel("Delta temperature in [K]")
    plt.title("Axial delta temp profile of fuel pellet (inner and outer) btw COLD and HOT")
    plt.legend()
    plt.grid()
    plt.show()

    # plot difference btw cladding hot and cold...
    plt.figure()
    plt.plot(xx,yy_cold_temp[:,1] - yy_hot_temp[:,1], label='Cladding external',color='red')
    plt.plot(xx,yy_cold_temp[:,2] - yy_hot_temp[:,2], label='Cladding internal',color='orange')
    plt.xlabel("Position in [m]")
    plt.ylabel("Delta temperature in [K]")
    plt.title("Axial delta temp profile of cladding (inner and outer) btw COLD and HOT")
    plt.legend()
    plt.grid()
    plt.show()


    #### ***************** PLOT TEMPERATURES (RADIAL) ******************* ####

    plt.figure()
    plt.plot(rr,rr_temp_fuel_radial[int( len(xx)/2 ),:], label="Temp profile")
    plt.plot(rr,np.ones(len(rr))*fuel_temp_max_suggested, label="Max suggested fuel temp", color='black', linestyle='--')
    plt.xlabel("Position in [m]")
    plt.ylabel("Temperature in [K]")
    plt.title("Middle position of pin, fuel pellet temperature profile")
    plt.grid()
    plt.show()

    # plot 3d temp radial profile of fuel
    x,y = np.meshgrid(xx,rr)
    fig_1 = plt.figure()
    ax = fig_1.add_subplot(111, projection='3d')
    ax.plot_surface(x,y,rr_temp_fuel_radial,cmap='viridis')
    ax.set_ylabel('Position along the pin [m]')
    ax.set_xlabel('Radius [m]')
    ax.set_zlabel('Temperature [K]')
    plt.show()


    #### ***************** OTHER PROPERTIES (AXIAL) ******************* ####
    # GAP
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

    # velocity of coolant along pin
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
    return None

plotting()