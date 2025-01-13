"""
*** HOT GEOMETRY THICKNESS AND MELTING MARGIN CALCULATION ***

"""

import matplotlib.pyplot as plt
from functions.thermal_functions import *
import os


burnup = 0
resol = 10
a = 0.455e-3
b = 0.535e-3
loadExisting = True


#### ******************************************* DOMAIN DISCRETIZATION ******************************************** ####
clad_range = np.arange(a,b,0.005e-3) # going to compute staff acc. to this clad thick. values
xx = np.linspace(pin_bottom_pos,pin_top_pos,resol)
yy_power_linear = np.zeros_like(xx)
yy_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
yy_hot_temp = np.zeros([len(xx),5]) # the same
yy_gap = np.zeros([len(xx),1])
yy_properties = np.zeros([len(xx),11])
clad_th = np.zeros_like(clad_range)
temp_max = np.zeros_like(clad_range)
margin_fuel = np.zeros_like(clad_range)
gap_min = np.zeros_like(clad_range)
gap_nom = np.zeros_like(clad_range)



#### **************************************************** CALCS *************************************************** ####
k = 0 # used in for cycle below

if loadExisting:
    gap_nom = np.load(os.path.join("main_thickCalc_numpy_saves", "gap_nom.npy"))
    clad_th = np.load(os.path.join("main_thickCalc_numpy_saves", "clad_th.npy"))
    temp_max = np.load(os.path.join("main_thickCalc_numpy_saves", "temp_max.npy"))
    margin_fuel = np.load(os.path.join("main_thickCalc_numpy_saves", "margin_fuel.npy"))
    gap_min = np.load(os.path.join("main_thickCalc_numpy_saves", "gap_min.npy"))
else:
    for clad_thick_0 in clad_range: # considering a value of cladding thickness from the array above...
        # COMPUTING
        for i in range(0,len(xx)): # working at z pos...
            yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i], yy_properties[i,:],clad_diam_out,fuel_diam_outer,_ = hot_geometry_general(
                xx[i], clad_d_outer, fuel_d_outer, clad_thick_0, burnup, print_status=False)

        gap_nom[k] = ( clad_d_outer/2 - clad_thick_0 - fuel_d_outer/2 )*1e6
        clad_th[k] = clad_thick_0 * 1000
        temp_max[k] = np.max(yy_hot_temp)
        margin_fuel[k] = fuel_temp_melting(burnup=burnup) - temp_max[k] # hp T melting lower
        gap_min[k] = np.min(yy_gap)*1e6

        print(f"Cladding thickness: {np.round(clad_th[k], 3)} mm")
        print(f"Max fuel inner temp: {np.round(temp_max[k], 2)} K")
        print(f"Margin fuel to melting: {np.round(margin_fuel[k], 2)} K")
        print(f"Gap: {np.round(gap_min[k], 3)} um")
        print(f"Gap nom: {gap_nom[k]} um")
        print(f"% gap: {100 * gap_min[k] / gap_nom[k]} %\n")

        k += 1


    np.save(os.path.join("main_thickCalc_numpy_saves", "gap_nom.npy"), gap_nom)
    np.save(os.path.join("main_thickCalc_numpy_saves", "clad_th.npy"), clad_th)
    np.save(os.path.join("main_thickCalc_numpy_saves", "temp_max.npy"), temp_max)
    np.save(os.path.join("main_thickCalc_numpy_saves", "margin_fuel.npy"), margin_fuel)
    np.save(os.path.join("main_thickCalc_numpy_saves", "gap_min.npy"), gap_min)



#### ************************************************* PLOTTING *************************************************** ####
plt.figure(1,figsize=(16, 9))

plt.plot(clad_th,temp_max,label='Max fuel temp')
plt.plot(xx, np.ones(len(xx)) * fuel_temp_melting(burnup=burnup), color='black', linestyle='--', label = 'Melting temp (CONS)') # tmelt provvisoria
plt.plot(xx, np.ones(len(xx)) * fuel_temp_max_suggested, color='black', linestyle='--', label = 'Max suggested temp') # tmelt provvisoria

plt.title("Max fuel temp acc. to clad thickness @ 0 GWd/ton")
plt.ylabel("K")
plt.xlabel("Cladding thickness value in [mm]")
plt.xlim([a*1000-0.005,b*1000+0.005])
plt.grid()
plt.savefig(os.path.join("figures","maxFuelTemp_vs_thickness_0.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



plt.figure(2,figsize=(16, 9))

plt.plot(clad_th,margin_fuel)

plt.title("Margin fuel temp to melting @ 0 GWd/ton")
plt.ylabel("Delta temp in K")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.savefig(os.path.join("figures","MarginFuelTempToMelt_0.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()



plt.figure(3,figsize=(16, 9))

plt.plot(clad_th,gap_min)

plt.title("Gap min acc. to cladding thickness @ 0 GWd/ton")
plt.ylabel("Gap [mm]")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.savefig(os.path.join("figures","GapMin_vs_Thickness_0.png"),dpi=300, bbox_inches='tight')
plt.show()
plt.close()