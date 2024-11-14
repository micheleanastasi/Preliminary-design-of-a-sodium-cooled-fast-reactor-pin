"""
CALLABLE FUNCTIONS
"""

from thermal_functions import *

def cold_geometry_temperature_computing(z):
    yy_temp_coolant = temp_coolant(z)
    yy_htc_local, yy_adim_number_cool, yy_cool_local_prop = heat_transfer_coefficient(yy_temp_coolant, clad_d_outer)
    yy_temp_clad_out = temp_cladding_outer(z, clad_d_outer)
    yy_temp_clad_in = temp_cladding_inner(z, clad_d_outer, clad_thickness_0)
    yy_temp_fuel_out, _ = temp_fuel_outer(z, clad_d_outer, fuel_d_outer, clad_thickness_0)
    yy_temp_fuel_in = temp_fuel_inner(z, clad_d_outer, fuel_d_outer, clad_thickness_0)

    output = [yy_temp_coolant,yy_temp_clad_out,yy_temp_clad_in,yy_temp_fuel_out,yy_temp_fuel_in]
    return output

