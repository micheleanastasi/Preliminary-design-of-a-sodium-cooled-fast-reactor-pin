"""
Useful script to demonstre negligible change of diameter after restructutring
"""
from functions.thermal_functions import *
import matplotlib.pyplot as plt

#### ********************************************* DOMAIN DISCRETIZATION ****************************************** ####
zz = np.linspace(pin_bottom_pos, pin_top_pos, 10)
rr = np.linspace(0,fuel_d_outer/2,15)
zz_power_linear = np.zeros_like(zz)
zz_cold_temp = np.zeros([len(zz), 5]) # coolant, clad out, clad in, fuel out, fuel in
zz_temp = np.zeros([len(zz), 5]) # the same
zz_r_clmn = np.zeros_like(zz)
zz_r_void = np.zeros_like(zz)

burnup=0

#### **************************************************** CALCS *************************************************** ####

# r clmn and r void
test = 0

temp_matrix = np.zeros([len(zz),len(rr)],dtype=float)


for i in range(0, len(zz)): # Z axis
    _, zz_temp[i, :], _, _, clad_d_out_hot, fuel_d_out_hot,_ = hot_geometry_general(zz[i],clad_d_outer,fuel_d_outer,clad_thickness_0,burnup,print_status=False)
    diam, old_diam, zz_r_clmn[i], zz_r_void[i], t_void, t_old = fuel_restructuring(zz[i], zz_temp[i, 3], zz_temp[i, 4], fuel_d_out_hot,clad_d_out_hot,burnup)

    print(f"Temp void: {np.round(t_void,2)} K - Temp old: {np.round(t_old,2)} K - DiamFuel {np.round(diam*1000,4)} mm - DiamOld {np.round(old_diam*1000,4)} mm")





plt.figure(1)
plt.plot(zz, zz_r_clmn)
plt.plot(zz,zz_r_void)
plt.grid()
plt.show()

plt.figure(2)
plt.plot(rr,temp_matrix[4,:])
#plt.plot(rr,temp_before[4,:])
plt.plot(rr,np.ones_like(rr)*fuel_temp_clmn,linestyle='--')
plt.grid()
plt.show()

#test = temp_before[4,:] - temp_matrix[4,:]
plt.plot(rr,test)
plt.grid()
plt.show()
