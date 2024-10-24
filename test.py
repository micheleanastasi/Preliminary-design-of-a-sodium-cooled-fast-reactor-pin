import numpy as np
from general_functions import *
from material_properties import *
from design_specifications import *
from thermal_functions import *

h = 130000
z = 0.85
clad_thick = 80e-6
temp_ci_guess = (650 + 273.15)
temp_cool = temp_coolant(z)
temp_co = temp_cladding_outer(temp_cool, power_lin_distribution(z), h)

print(clad_thermal_cond.subs(temp,873.15))


eqz_1 = temp - temp_co
print(eqz_1)
eqz_2 = power_lin_distribution(z) * clad_thick / (pi * (clad_d_outer - 2 * clad_thick) * clad_thermal_cond)
res = eqz_1 - eqz_2

out = iterative_solver(res,temp_ci_guess)
print(f"\n")
print(temp_co - 273.15)
print(out-273.15)