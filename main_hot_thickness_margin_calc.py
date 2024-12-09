"""
*** HOT GEOMETRY THICKNESS AND MELTING MARGIN CALCULATION ***

Main script to get data (only plots currently) in order to estimate an optimal cladding thickness according to gap and
margin to fuel melting (temp)
TIME: HOT GEO; HENCE FIRST MOMENTS (NO RESTR, REDISTR, BURN_UP)

NOTE: still missing phenomena such as restructuring, burn-up and redistribution (the latter has minor impact though)
    Nonetheless, up to this point everything was taken into account in a conservative approach, although some aspect should
    be evaluated accordingly (i.e. how much is Pu redistr. negligible?)
    But why conservative? Burn up -> inf, then lowering the values of T_fuel_melting, k_thermal, etc at minimum... moreover not still
    considering restructuring!
    THIS IS THE FIRST STEP

SO:
- start: only hot geometry, no burn up (we are HERE but cons.hp: burn up --> inf, to have K_fuel and Tm as low as possible, CONSERVATIVE!)
- then: hot geo + restructuring + Pu redistribution
- end: hot geo + "" + "" + burn up
properties changing: k_fuel, Tm_fuel, ...

NB CONTACT HEAT TRANSFER BTW FUEL AND CLADDING NOT CONSIDERED AS THERE MUST BE ALWAYS SPACE AT THIS STAGE!
"""

import numpy as np
import matplotlib.pyplot as plt
from thermal_functions import *


#### ******************************************* DOMAIN DISCRETIZATION ******************************************** ####
resol = 10
xx = np.linspace(pin_bottom_pos,pin_top_pos,resol)
rr = np.linspace(0,fuel_d_outer/2,resol)

clad_range = np.arange(0.485e-3,0.535e-3,0.005e-3) # going to compute staff acc. to this clad thick. values

# arrays initialization
yy_power_linear = np.zeros_like(xx)
yy_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
yy_hot_temp = np.zeros([len(xx),5]) # the same
yy_gap = np.zeros([len(xx),1])
yy_properties = np.zeros([len(xx),11])
rr_temp_fuel_radial = np.zeros((len(xx),len(rr)))
clad_th = np.zeros_like(clad_range)
temp_max = np.zeros_like(clad_range)
margin_fuel = np.zeros_like(clad_range)
gap_min = np.zeros_like(clad_range)
gap_nom = np.zeros_like(clad_range)

burnup = 0


#### **************************************************** CALCS *************************************************** ####
k = 0 # used in for cycle below
for clad_thick_0 in clad_range: # considering a value of cladding thickness from the array above...
    test = 0 # used to keep tracking the computing process
    gap_nom[k] = ( clad_d_outer/2 - clad_thick_0 - fuel_d_outer/2 )*1e6
    for i in range(0,len(xx)): # working at z pos...
        yy_power_linear[i] = power_lin_distribution(xx[i])  # linear power calc at z pos
        yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i], yy_properties[i,:],clad_diam_out,fuel_diam_outer = hot_geometry_general(
            xx[i], clad_d_outer, fuel_d_outer, clad_thick_0, burnup, print_status=False)
        #for j in range(0,len(rr)):
        #    rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j], xx[i], clad_diam_out, fuel_diam_outer, clad_thick_0)
        #    test += 1
    clad_th[k] = clad_thick_0 * 1000
    temp_max[k] = np.max(yy_hot_temp)
    margin_fuel[k] = fuel_temp_melting - temp_max[k] # hp T melting lower
    gap_min[k] = np.min(yy_gap)*1e6
    print(f"Cladding thickness: {np.round(clad_th[k],3)} mm")
    print(f"Max fuel inner temp: {np.round(temp_max[k],2)} K")
    print(f"Margin fuel to melting: {np.round(margin_fuel[k],2)} K")
    print(f"Gap: {np.round(gap_min[k],3)} um")
    print(f"Gap nom: {gap_nom[k]} um")
    print(f"% gap: {100*gap_min[k]/gap_nom[k]} %\n")
    k += 1



#### ************************************************* PLOTTING *************************************************** ####
plt.figure(1)
plt.plot(clad_th,temp_max,label='Max fuel temp')
plt.xlim([0.47,0.55])
plt.plot(xx, np.ones(len(xx)) * 2964.92, color='black', linestyle='--', label = 'Melting temp (CONS)') # tmelt provvisoria
plt.plot(xx, np.ones(len(xx)) * fuel_temp_max_suggested, color='black', linestyle='--', label = 'Max suggested temp') # tmelt provvisoria
plt.title("Max fuel temp acc. to clad thickness")
plt.ylabel("K")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()

plt.figure(2)
plt.plot(clad_th,margin_fuel)
plt.title("Margin fuel temp to melting")
plt.ylabel("Delta temp in K")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()

plt.figure(3)
plt.plot(clad_th,gap_min)
plt.title("Gap min acc. to cladding thickness")
plt.ylabel("Gap [mm]")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()