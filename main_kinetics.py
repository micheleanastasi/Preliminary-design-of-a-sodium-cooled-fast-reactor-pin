
import numpy as np
from thermal_functions import *
import matplotlib.pyplot as plt

bb = np.arange(0,52,4)
xx = np.linspace(pin_bottom_pos, pin_top_pos, 5)

yy_hot_temp = np.zeros([5])
yy_gaps = np.zeros_like(bb)
plot_temp = np.zeros_like(bb)

for j in range(0,len(bb)):
    _,yy_hot_temp[:], yy_gaps[j],_,_,_ = hot_geometry_general(0.45, clad_d_outer, fuel_d_outer, clad_thickness_0,bb[j],print_status=True)
    plot_temp[j] = yy_hot_temp[4]
    print(f"Burnup: {bb[j]} GWd/ton")

plt.figure()
plt.plot(bb,plot_temp)
plt.grid()
plt.show()

plt.figure()
plt.plot(bb,yy_gaps)
plt.grid()
plt.show()
"""
yy_cold_temp[i,:], yy_hot_temp[i,:], yy_gap[i], yy_properties[i,:], clad_diam_out[i], fuel_diam_outer[i] = hot_geometry_general(
        xx[i], clad_d_outer, fuel_d_outer, clad_thickness_0,bup=burnup)

"""