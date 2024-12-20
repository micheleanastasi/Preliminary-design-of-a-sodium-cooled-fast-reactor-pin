"""
GENERAL MAIN SCRIPT

KINETICS:
- start: only hot geometry, no burn up (we are HERE but cons.hp: burn up --> inf)
- then: hot geo + restructuring
- end: hot geo + "" + burn up effects (then swelling)

PROPERTIES changing w.r.t burnup: k_fuel, Tm_fuel, k_gas
"""

# IMPORTING
import os
import matplotlib.pyplot as plt
import numpy as np

from functions.thermal_functions import *

## ADJUSTABLE PARAMETERS
burnup = (0,1,52,104) # GWd/ton
res = 40
extra_pin_len = 0.75 # m - little diameter expansion then (whereas length exp neglected!) (HP CONS) !
# thickness from general_properties.py : 0.48 mm !
# initial gap size: 85 um !

loadExisting = False



## DOMAIN
xx = np.linspace(pin_bottom_pos, pin_top_pos, res)

## OUTPUTS
yy_power_linear = np.zeros_like(xx) # constant all the time!
yy_cold_temp = np.zeros([len(xx),5,len(burnup)]) # coolant, clad out, clad in, fuel out, fuel in

yy_hot_temp = np.zeros([len(xx),5,len(burnup)]) # the same
yy_gap = np.zeros([len(xx),len(burnup)])
yy_properties = np.zeros([len(xx),11,len(burnup)])
clad_diam_out = np.zeros([len(xx),len(burnup)])
fuel_diam_outer = np.zeros([len(xx),len(burnup)])

yy_r_clmn = np.zeros([len(xx),len(burnup)])
yy_r_void = np.zeros([len(xx),len(burnup)])
yy_temp_void = np.zeros([len(xx),len(burnup)])

vol_hot = np.zeros(len(burnup))
pressure = np.zeros(len(burnup))
extra_pin = np.zeros(len(burnup))



if loadExisting:
    yy_hot_temp = np.load(os.path.join("main_numpy_saves", "yy_hot_temp.npy"))
    yy_cold_temp = np.load(os.path.join("main_numpy_saves", "yy_cold_temp.npy"))
    yy_gap = np.load(os.path.join("main_numpy_saves", "yy_gap.npy"))
    yy_properties = np.load(os.path.join("main_numpy_saves", "yy_properties.npy"))
    clad_diam_out = np.load(os.path.join("main_numpy_saves", "clad_diam_out.npy"))
    fuel_diam_outer = np.load(os.path.join("main_numpy_saves", "fuel_diam_outer.npy"))
    vol_hot = np.load(os.path.join("main_numpy_saves", "vol_hot.npy"))
    pressure = np.load(os.path.join("main_numpy_saves", "pressure.npy"))
    extra_pin = np.load(os.path.join("main_numpy_saves", "extra_pin.npy"))
    yy_r_clmn = np.load(os.path.join("main_numpy_saves", "yy_r_clmn.npy"))
    yy_r_void = np.load(os.path.join("main_numpy_saves", "yy_r_void.npy"))
    yy_temp_void = np.load(os.path.join("main_numpy_saves", "yy_temp_void.npy"))
else:
    ## COMPUTING AND THEN SAVING IN .NPY FORMAT
    for j in range(0,len(burnup)):
        print(f"\n\n\n\n*************** BURNUP = {burnup[j]} GWd/ton *****************************************************")
        for i in range(0,len(xx)): # Z axis
            yy_power_linear[i] = power_lin_distribution(xx[i])

            # COMPUTING TEMPS, GAPS, NEW CLAD AND FUEL DIAMETERS
            yy_cold_temp[i,:,j], yy_hot_temp[i, :,j], yy_gap[i,j], yy_properties[i,:,j], clad_diam_out[i,j], fuel_diam_outer[i,j] = hot_geometry_general(
                xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0, burnup[j])

            # RESTRUCTURING
            _, _, yy_r_clmn[i,j], yy_r_void[i,j], yy_temp_void[i,j], _ = fuel_restructuring( xx[i], yy_hot_temp[i,3,j], yy_hot_temp[i,4,j], fuel_diam_outer[i,j],
                                                                                             clad_diam_out[i,j], burnup[j] )

        ## PRESSURE COMPUTING
        vol_hot[j] = gap_vol_hot(fuel_diam_outer[:,j], clad_diam_out[:,j]) # arrays as input
        # computed considering expanded clad diam (inner) as hot geometry
        vol_plenum = extra_pin_len * 0.25 * pi * clad_d_inner ** 2
        pressure[j], _ = pressure_gap_calc( vol_hot[j], yy_hot_temp[:,2:4,j], burnup[j], vol_plenum,
                                            plenum_clad_d_in=clad_d_inner, temp_plenum=yy_hot_temp[0, 0,j],print_stuff=True )

    # SAVING...
    np.save(os.path.join("main_numpy_saves","yy_hot_temp.npy"),yy_hot_temp)
    np.save(os.path.join("main_numpy_saves","yy_cold_temp.npy"),yy_cold_temp)
    np.save(os.path.join("main_numpy_saves","yy_gap.npy"),yy_gap)
    np.save(os.path.join("main_numpy_saves","yy_properties.npy"),yy_properties)
    np.save(os.path.join("main_numpy_saves","clad_diam_out.npy"),clad_diam_out)
    np.save(os.path.join("main_numpy_saves","fuel_diam_outer.npy"),fuel_diam_outer)
    np.save(os.path.join("main_numpy_saves","vol_hot.npy"),vol_hot)
    np.save(os.path.join("main_numpy_saves","pressure.npy"),pressure)
    np.save(os.path.join("main_numpy_saves","extra_pin.npy"),extra_pin)
    np.save(os.path.join("main_numpy_saves","yy_r_clmn.npy"),yy_r_clmn)
    np.save(os.path.join("main_numpy_saves","yy_r_void.npy"),yy_r_void)
    np.save(os.path.join("main_numpy_saves","yy_temp_void.npy"),yy_temp_void)



## ******************************************** PLOT **************************************************************** ##
## Axial temp profile of fuel pellet (inner and outer) ##
plt.figure(1,figsize=(16, 9))

plt.plot(xx,yy_hot_temp[:,3,0], label='Fuel external @ 0 GWd/ton',color='blue',linestyle='-')
plt.plot(xx,yy_hot_temp[:,3,1], label='Fuel external @ 1 GWd/ton',color='blue',linestyle='--')
plt.plot(xx,yy_hot_temp[:,3,2], label='Fuel external @ 52 GWd/ton',color='blue',linestyle=':')

plt.plot(xx,yy_hot_temp[:,4,0], label='Fuel external @ 0 GWd/ton',color='red',linestyle='-')
plt.plot(xx,yy_hot_temp[:,4,1], label='Fuel external @ 1 GWd/ton',color='red',linestyle='--')
plt.plot(xx,yy_hot_temp[:,4,2], label='Fuel external @ 52 GWd/ton',color='red',linestyle=':')

plt.plot(xx,np.ones(len(xx))*fuel_temp_max_suggested, label='Max suggested fuel temp', color='gray', linestyle='--')
plt.plot(xx,np.ones(len(xx))*fuel_temp_melting(burnup=burnup[0]), label='Melting point of fuel @ 0 GWd/ton)', color='black', linestyle='-')
plt.plot(xx,np.ones(len(xx))*fuel_temp_melting(burnup=burnup[1]), label='Melting point of fuel @ 1 GWd/ton)', color='black', linestyle='--')
plt.plot(xx,np.ones(len(xx))*fuel_temp_melting(burnup=burnup[2]), label='Melting point of fuel @ 52 GWd/ton)', color='black', linestyle=':')

plt.xlabel("Position along the pin in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temp profile of fuel pellet (inner and outer)")
plt.legend(loc='best')
plt.grid()
plt.savefig(os.path.join("saved_figures","fuel_pellet_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## plot coolant, cladding in and out temp for hot and cold geometries
plt.figure(2,figsize=(16, 9))

plt.plot(xx,yy_cold_temp[:,1,0], label='COLD Cladding external',color='red', linestyle='--')
plt.plot(xx,yy_cold_temp[:,2,0], label='COLD Cladding internal',color='orange', linestyle='--')
plt.plot(xx,yy_hot_temp[:,0,0], label='Coolant',color='blue')
plt.plot(xx,yy_hot_temp[:,1,0], label='HOT Cladding external',color='red')
plt.plot(xx,yy_hot_temp[:,2,0], label='HOT Cladding internal',color='orange')
plt.plot(xx,np.ones(len(xx))*clad_temp_max, label='Max suggested cladding temp', color='black', linestyle='--')

plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temp profile of coolant and cladding (inner and outer) @ 0 GWd/ton")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","coolClad_0.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## Axial temperature difference between cladding inner e cladding outer) w.r.t burn up ##
plt.figure(3,figsize=(16, 9))

delta_temp_cl_0 = yy_hot_temp[:,2,0] - yy_hot_temp[:,1,0]
delta_temp_cl_1 = yy_hot_temp[:,2,1] - yy_hot_temp[:,1,1]
delta_temp_cl_52 = yy_hot_temp[:,2,2] - yy_hot_temp[:,1,2]
plt.plot(xx,delta_temp_cl_0, label='Difference @ 0 GWd/ton')
plt.plot(xx,delta_temp_cl_1, label='Difference @ 1 GWd/ton')
plt.plot(xx,delta_temp_cl_52, label='Difference @ 52 GWd/ton')

plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temperature difference between cladding inner e cladding outer) w.r.t burn up")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","deltaCladIn_vs_CladOut_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## Axial temperature difference between fuel outer e cladding inner w.r.t burn up ##
plt.figure(4,figsize=(16, 9))

delta_temp_cl_0 = yy_hot_temp[:,3,0] - yy_hot_temp[:,2,0]
delta_temp_cl_1 = yy_hot_temp[:,3,1] - yy_hot_temp[:,2,1]
delta_temp_cl_52 = yy_hot_temp[:,3,2] - yy_hot_temp[:,2,2]
plt.plot(xx,delta_temp_cl_0, label='Difference @ 0 GWd/ton')
plt.plot(xx,delta_temp_cl_1, label='Difference @ 1 GWd/ton')
plt.plot(xx,delta_temp_cl_52, label='Difference @ 52 GWd/ton')

plt.xlabel("Position in [m]")
plt.ylabel("Temperature in [K]")
plt.title("Axial temperature difference between fuel outer e cladding inner w.r.t burn up")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","deltaFuelOuter_vs_CladInner_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## Fuel diameter variation w.r.t burnup ##
plt.figure(5,figsize=(16, 9))

plt.plot(xx,fuel_diam_outer[:,0]*1e3, label='Fuel diameter # 0 GWd/ton')
plt.plot(xx,fuel_diam_outer[:,1]*1e3, label='Fuel diameter # 1 GWd/ton')
plt.plot(xx,fuel_diam_outer[:,2]*1e3, label='Fuel diameter # 52 GWd/ton')

plt.plot(xx,np.ones(len(xx))*fuel_d_outer*1e3, label='Initial fuel diameter', color='black', linestyle='--')

plt.xlabel("Position along the pin in [m]")
plt.ylabel("Diameter in [mm]")
plt.title("Fuel diameter variation w.r.t burnup")
plt.legend(loc='best')
plt.grid()
plt.savefig(os.path.join("saved_figures","fuel_diameter_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## Cladding diameter variation w.r.t burnup ##
plt.figure(6,figsize=(16, 9))

plt.plot(xx,clad_diam_out[:,0]*1e3, label='Cladding diameter # 0 GWd/ton')
plt.plot(xx,clad_diam_out[:,1]*1e3, label='Cladding diameter # 1 GWd/ton')
plt.plot(xx,clad_diam_out[:,2]*1e3, label='Cladding diameter # 52 GWd/ton')

plt.plot(xx,np.ones(len(xx))*clad_d_outer*1e3, label='Initial fuel diameter', color='black', linestyle='--')

plt.xlabel("Position along the pin in [m]")
plt.ylabel("Diameter in [mm]")
plt.title("Cladding diameter variation w.r.t burnup")
plt.legend(loc='best')
plt.grid()
plt.savefig(os.path.join("saved_figures","clad_diameter_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## gap size variation w.r.t burnup ##
plt.figure(7,figsize=(16, 9))

plt.plot(xx,yy_gap[:,0]*1e6, label='Gap size # 0 GWd/ton')
plt.plot(xx,yy_gap[:,1]*1e6, label='Gap size # 1 GWd/ton')
plt.plot(xx,yy_gap[:,2]*1e6, label='Gap size # 52 GWd/ton')

plt.plot(xx,np.ones(len(xx))*initial_delta_gap*1e6, label='Initial delta gap', color='black', linestyle='--')

plt.xlabel("Position along the pin in [m]")
plt.ylabel("Gap size in [um]")
plt.title("Gap size variation w.r.t burnup")
plt.legend(loc='best')
plt.grid()
plt.savefig(os.path.join("saved_figures","gapSize_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## radius of void and columnar
plt.figure(8,figsize=(16, 9))

plt.plot(xx,yy_r_clmn[:,0]*1e3,label="Columnar region radius")
plt.plot(xx,yy_r_void[:,0]*1e3,label="Void region radius")
plt.plot(xx,fuel_diam_outer[:,0]*1e3,label="Fuel pellet radius")

plt.xlabel("Position along the pin in [m]")
plt.ylabel("Radius size in [mm]")
plt.title("Restructuring effects: radius of regions")
plt.legend(loc='best')
plt.grid()
plt.savefig(os.path.join("saved_figures","rVoid_rClmn_0.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## pressure ##
plt.figure(9,figsize=(16, 9))

plt.plot(burnup[0:3],pressure[0:3]*1e-6)

plt.xlabel("Burn up in [GWd/ton]")
plt.ylabel("Pressure in [MPa]")
plt.title("Pressure variation w.r.t burnup")
plt.grid()
plt.savefig(os.path.join("saved_figures","pressure_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()


## velocity
max_vel = 8 # m/s
plt.figure(10,figsize=(16, 9))

plt.plot(xx,yy_properties[:,5], label='Velocity')
plt.plot(xx,np.ones_like(xx)*max_vel, label='Maximum velocity allowed', linestyle='--', color='black')

plt.xlabel("Position [m]")
plt.ylabel("Velocity [m/s]")
plt.title("Average velocity (over the z section)")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","coolVel_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## margin
plt.figure(11,figsize=(16, 9))

margin_0 = fuel_temp_melting(burnup=burnup[0]) - yy_hot_temp[:,4,0]
margin_1 = fuel_temp_melting(burnup=burnup[1]) - yy_hot_temp[:,4,1]
margin_52 = fuel_temp_melting(burnup=burnup[2]) - yy_hot_temp[:,4,2]

plt.plot(xx,margin_0, label='Margin @ 0 GWd/ton')
plt.plot(xx,margin_1, label='Margin @ 1 GWd/ton')
plt.plot(xx,margin_52, label='Margin @ 52 GWd/ton')

plt.xlabel("Position [m]")
plt.ylabel("Temperature in [K]")
plt.title("Margin variation w.r.t. Burn Up along the pin")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","marginToFuelMelt_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## cladding temperature
plt.figure(12,figsize=(16, 9))

plt.plot(xx,yy_hot_temp[:,1,0], label='Cladding external temp @ 0 GWd/ton',color='blue',linestyle='-')
plt.plot(xx,yy_hot_temp[:,1,1], label='Cladding external temp @ 1 GWd/ton',color='blue',linestyle='--')
plt.plot(xx,yy_hot_temp[:,1,2], label='Cladding external temp @ 52 GWd/ton',color='blue',linestyle=':')

plt.plot(xx,yy_hot_temp[:,2,0], label='Cladding internal temp @ 0 GWd/ton',color='red',linestyle='-')
plt.plot(xx,yy_hot_temp[:,2,1], label='Cladding internal temp @ 1 GWd/ton',color='red',linestyle='--')
plt.plot(xx,yy_hot_temp[:,2,2], label='Cladding internal temp @ 52 GWd/ton',color='red',linestyle=':')

plt.xlabel("Position [m]")
plt.ylabel("Temperature in [K]")
plt.title("Cladding temperature variation w.r.t. Burn Up along the pin")
plt.legend()
plt.grid()
plt.savefig(os.path.join("saved_figures","claddingTemp_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



## pressure (extra) - calcolo pressione in funzione di extra volume (in termini di lunghezza)
vol_extra = np.arange(0.5e-6,55e-6,0.5e-6)
res_y = np.zeros([len(vol_extra),len(burnup)])
res_x = np.zeros([len(vol_extra),len(burnup)])

for j in range(0,len(burnup[0:3])):
    c = 0
    for i in vol_extra:
        # calcolo considering initial coolant temperature (HP SEMPL) -- yy_hot_temp[0,0,...]
        res_y[c,j], res_x[c,j] = pressure_gap_calc(vol_hot[j], yy_hot_temp[:,2:4,0], burnup[j], plenum_vol=i, plenum_clad_d_in=clad_d_inner,
                                               temp_plenum=yy_hot_temp[0, 0,j],print_stuff=False)
        c += 1


plt.figure(13,figsize=(16, 9))

plt.plot(res_x[:,0]*1000,res_y[:,0]/1e6,label='Pressure @ 0 GWd/ton')
plt.plot(res_x[:,0]*1000,res_y[:,1]/1e6,label='Pressure @ 1 GWd/ton')
plt.plot(res_x[:,0]*1000,res_y[:,2]/1e6,label='Pressure @ 52 GWd/ton')
plt.plot(res_x[:,0]*1000,np.ones_like(res_x[:,0])*5,label='Maximum suggested pressure (5 MPa)',color='black',linestyle='--')

plt.plot(res)
plt.yscale("log")

plt.xlabel("Extra pin length in [mm]")
plt.ylabel("Total pressure in cladding in [MPa]")
plt.title("Pressure variation w.r.t burnup and extra length (volume) added")
plt.grid()
plt.legend()
plt.savefig(os.path.join("saved_figures","press_vs_extraLength_0_1_52.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()