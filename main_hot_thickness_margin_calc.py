"""
*** HOT GEOMETRY THICKNESS AND MELTING MARGIN CALCULATION ***

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

NB CONTACT HEAT TRANSFER BTW FUEL AND CLADDING NOT CONSIDERED
"""

import numpy as np
import matplotlib.pyplot as plt
from thermal_functions import *


#### ***************** DOMAIN DISCRETIZATION ********************** ####
resol = 15
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

    #print(f"Hot geo completed at {np.round(100*z/0.85,2)}% (Position: {np.round(z,2)} m) - Temp fuel inner: HOT:{np.round(temp_array[4],2)}, COLD:{np.round(old[4],2)} - Gap:{np.round(delta_gap*1000,6)} mm")
    other = np.array( list([yy_htc_loc]) + list(yy_adim_num_cool) + list(yy_cool_loc_prop) )
    return old,temp_array,delta_gap,other, clad_d_out_0, fuel_d_out_0



#### *********************** CALCS *********************** ####

yy_power_linear = np.zeros_like(xx)
yy_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
yy_hot_temp = np.zeros([len(xx),5]) # the same
yy_gap = np.zeros([len(xx),1])
yy_properties = np.zeros([len(xx),11])
rr_temp_fuel_radial = np.zeros((len(xx),len(rr)))

clad_range = np.arange(0.505e-3,0.570e-3,0.005e-3)
clad_th = np.zeros_like(clad_range)
temp_max = np.zeros_like(clad_range)
margin_fuel = np.zeros_like(clad_range)
gap_min = np.zeros_like(clad_range)
k = 0

for clad_thickness_0 in clad_range:
    test = 0
    for i in range(0,len(xx)):
        yy_power_linear[i] = power_lin_distribution(xx[i])
        yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i],yy_properties[i,:],clad_diam_out,fuel_diam_outer = hot_geometry_iteration(xx[i],clad_d_outer,fuel_d_outer,clad_thickness_0)
        for j in range(0,len(rr)):
            rr_temp_fuel_radial[i,j] = temp_fuel_inner_radial(rr[j],xx[i],clad_diam_out,fuel_diam_outer,clad_thickness_0)
            test += 1
            #print(f"Calcs completed at {np.round(100*test/(len(xx)*len(rr)),3)} %")
    clad_th[k] = clad_thickness_0*1000
    temp_max[k] = np.max(yy_hot_temp)
    margin_fuel[k] = fuel_temp_melting - temp_max[k] # hp T melting lower
    gap_min[k] = np.min(yy_gap)*1e6
    print(f"Cladding thickness: {np.round(clad_th[k],3)} mm")
    print(f"Max fuel inner temp: {np.round(temp_max[k],2)} K")
    print(f"Margin fuel to melting: {np.round(margin_fuel[k],2)} K")
    print(f"Gap: {np.round(gap_min[k],3)} um\n")
    k += 1

plt.figure()
plt.plot(clad_th,temp_max,label='Max fuel temp')
plt.xlim([0.50,0.57])
plt.plot(xx, np.ones(len(xx)) * 2964.92, color='black', linestyle='--', label = 'Melting temp (CONS)') # tmelt provvisoria
plt.plot(xx, np.ones(len(xx)) * fuel_temp_max_suggested, color='black', linestyle='--', label = 'Max suggested temp') # tmelt provvisoria
plt.title("Max fuel temp acc. to clad thickness")
plt.ylabel("K")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()

plt.figure()
plt.plot(clad_th,margin_fuel)
plt.title("Margin fuel temp to melting")
plt.ylabel("Delta temp in K")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()

plt.plot(clad_th,gap_min)
plt.title("Gap min acc. to cladding thickness")
plt.ylabel("Gap [mm]")
plt.xlabel("Cladding thickness value in [mm]")
plt.grid()
plt.show()