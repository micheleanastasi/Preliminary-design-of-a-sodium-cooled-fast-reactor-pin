
from thermal_functions import *


#### ********************************************* DOMAIN DISCRETIZATION ****************************************** ####
xx = np.linspace(pin_bottom_pos,pin_top_pos,20)
rr = np.linspace(0,fuel_d_outer/2,10)
zz_power_linear = np.zeros_like(xx)
zz_cold_temp = np.zeros([len(xx),5]) # coolant, clad out, clad in, fuel out, fuel in
zz_hot_temp = np.zeros([len(xx), 5]) # the same
zz_r_clmn = np.zeros_like(xx)
zz_r_void = np.zeros_like(xx)


#### **************************************************** CALCS *************************************************** ####

# r clmn and r void
test = 0
for i in range(0,len(xx)): # Z axis
    z = xx[i]
    zz_power_linear[i] = power_lin_distribution(xx[i])
    _, zz_hot_temp[i, :], _, _, clad_d_out_hot, fuel_d_out_hot = hot_geometry_iteration(xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)

    fuel_r_out_hot = fuel_d_out_hot/2
    zz_r_clmn[i] = get_R_from_temp(z,fuel_r_out_hot,fuel_temp_clmn,zz_hot_temp[i,4])
    zz_r_void[i] = radius_void_get( zz_r_clmn[i], poro_asf )