from numpy.ma.core import zeros_like

from thermal_functions import *
import matplotlib.pyplot as plt

tt = np.linspace(0, 100, 1)  # MWd/kg # RIVEDI
def max_temp_kinetics():
    xx = np.linspace(pin_bottom_pos, pin_top_pos, 20)


    out = np.zeros_like(tt)
    yy_power_linear = np.zeros_like(xx)
    yy_cold_temp = np.zeros([len(xx), 5])  # coolant, clad out, clad in, fuel out, fuel in
    yy_hot_temp = np.zeros([len(xx), 5])  # the same
    yy_gap = np.zeros([len(xx), 1])
    yy_properties = np.zeros([len(xx), 11])

    for t in range(0, len(tt)):
        for i in range(0, len(xx)):  # Z axis
            yy_power_linear[i] = power_lin_distribution(xx[i])
            yy_cold_temp[i, :], yy_hot_temp[i, :], yy_gap[i], yy_properties[i,:], clad_diam_out, fuel_diam_outer = hot_geometry_iteration(
                xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0)
        out[t] = np.max( yy_hot_temp[:,:], axis=0 )
    return out

temp = max_temp_kinetics()
plt.plot(tt,temp)