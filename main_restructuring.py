import numpy as np

from thermal_functions import *
import matplotlib.pyplot as plt

#### ********************************************* DOMAIN DISCRETIZATION ****************************************** ####
zz = np.linspace(pin_bottom_pos, pin_top_pos, 10)
rr = np.linspace(0,fuel_d_outer/2,15)
zz_power_linear = np.zeros_like(zz)
zz_cold_temp = np.zeros([len(zz), 5]) # coolant, clad out, clad in, fuel out, fuel in
zz_hot_temp = np.zeros([len(zz), 5]) # the same
zz_r_clmn = np.zeros_like(zz)
zz_r_void = np.zeros_like(zz)


#### **************************************************** CALCS *************************************************** ####

# r clmn and r void
test = 0

temp_matrix = np.zeros([len(zz),len(rr)])
temp_before = np.zeros([len(zz),len(rr)])

for i in range(0, len(zz)): # Z axis
    # get R columnar, R void for each point of discrete domain along z
    z = zz[i]
    zz_power_linear[i] = power_lin_distribution(zz[i])
    _, zz_hot_temp[i, :], _, _, clad_d_out_hot, fuel_d_out_hot = hot_geometry_iteration(zz[i], clad_d_outer, fuel_d_outer, clad_thickness_0)

    fuel_r_out_hot = fuel_d_out_hot/2
    zz_r_clmn[i] = get_R_from_temp(z,fuel_r_out_hot,fuel_temp_clmn,zz_hot_temp[i,4],zz_hot_temp[i,3])
    zz_r_void[i] = radius_void_get( zz_r_clmn[i], poro_asf )

    temp_out = zz_hot_temp[i, 3]  # temp @ radius fuel outer - it's the same in spite of restr.
    r_c = zz_r_clmn[i]
    r_v = zz_r_void[i]

    ## NOTE x, pu changing due to kinetics and redistribution, in order
    ## firstly, temp variation due to restructuring computing
    for j in range(len(rr)-1, -1, -1): # da len-1 to 0 index
        temp_before[i, j] = temp_fuel_inner_radial(rr[j], zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0)
        if rr[j] >= r_v:
            temp_matrix[i,j] = temp_fuel_inner_radial(rr[j], zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0,x=0,pu=0.29,po=0.12,r_c=r_c,r_v=r_v)
        else:
            temp_matrix[i, j] = temp_fuel_inner_radial(r_v, zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0, x=0, pu=0.29, po=0.12, r_c=r_c, r_v=r_v)
    ## ancora da implementare regioni asf e clmn, calcolo solo T max (void factor) radiale
"""
# WIP
        if r_c == 0 and r_v == 0: # zones where there isn't restr.
            temp_matrix[i,j] = temp_fuel_inner_radial(rr[j], zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0)
        else: # restr. computing
            if rr[j] >= r_c: # as fabricated
                temp_matrix[i,j] = temp_fuel_inner_radial(rr[j], zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0,x=0,pu=0.29,po=0.12)
                last_one = temp_matrix[i,j]
            elif r_c > rr[j] > r_v: # columnar grains region
                temp_matrix[i,j] = temp_fuel_inner_radial(rr[j], zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0,x=0,pu=0.29,po=0.05,r_c=r_c,r_v=r_v)
            elif rr[j] <= r_v: # void region
                temp_matrix[i,j] = temp_fuel_inner_radial(r_v, zz[i], clad_d_out_hot, fuel_d_out_hot, clad_thickness_0,x=0,pu=0.29,po=0.05,r_c=r_c,r_v=r_v)
"""


plt.figure(1)
plt.plot(zz, zz_r_clmn)
plt.plot(zz,zz_r_void)
plt.plot(zz,np.ones_like(zz)*fuel_d_out_hot/2)
plt.grid()
plt.show()

plt.figure(2)
plt.plot(rr,temp_matrix[4,:])
plt.plot(rr,temp_before[4,:])
plt.plot(rr,np.ones_like(rr)*fuel_temp_clmn,linestyle='--')
plt.grid()
plt.show()

test = temp_before[4,:] - temp_matrix[4,:]
plt.plot(rr,test)
plt.grid()
plt.show()