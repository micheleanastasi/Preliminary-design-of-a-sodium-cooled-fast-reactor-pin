from general_functions import *
from material_properties import *
from design_specifications import *

"""
#### NOTE ####
#
# c'è ancora da considerare burn up effects
# anche fission gases --> K_gas
#
"""


#### *********** TEMP PROFILE OF COOLANT ALONG Z AXIS ************ ####
#
# energy balance: f*int(Cp(T), Tin, T) - int(q', bottom_pin, z) / mass_flow = 0 ... Tout is our goal.
#
## occorre metodo che minimizzi la funzione res affinche si trovi una temperatura che rispetti il bilancio energetico, in quanto Cp stesso dipende da temperatura!!
def temp_coolant(z):
    """
    :param z: position choice
    :return: Temp of coolant at point T
    """
    T_0 = 800 + 273.15 # K
    eqz_2 = f_pitch*integral_power_lin_distr(z) / cool_mass_flow  # valore num seconda parte eqz !! NB F PITCH=1/2 !!
    eqz_1 = sy.integrate(cool_spec_heat,temp)
    res = eqz_1 - eqz_1.subs(temp,cool_temp_inlet) - eqz_2 # primo è integtale, secondo intehrale con sostituzione (quindi risolto) e terzo eqz_2

    T_coolant_at_z = equation_temp_solver(res, T_0)
    return T_coolant_at_z


## CORRELATION FUNCTIONS
#
# Calculated for each temperature of the coolant along z axis
# HP: pitch btw pins constant, no cross-section of flows!!
#
#* {NB: possibile "dubbio" su htc, occorre verifica h e Pr}
def heat_transfer_coeff_local(temperature,clad_diam_out):
    """
    Get parameters of coolant such as heat transfer coeff., adim. numbers and avg_velocity: point by point, according to temperature and clad diameter!
    :param clad_diam_out: in m - Varying parameter! (hot geometry)
    :param temperature: in [K]
    :return: htc, array with adim. numbers and array with some coolant/geometry properties calculated here
    """
    # parameters that are going to be useful below...
    temp_fahr = (temperature - 273.15) * 9 / 5 + 32 # conversione kelvin to fahrenheit
    density = cool_density__fahr.subs(temp_f,temp_fahr) # kg/m^3
    dyn_viscosity = cool_dyn_viscosity.subs(temp,temperature) # Pa*s
    spec_heat = cool_spec_heat.subs(temp,temperature) # J/Kg/K
    th_cond = cool_th_cond.subs(temp,temperature) # J/m/K

    net_area = pin_pitch**2 * sqrt(3)/4 - 0.5*( pi*clad_diam_out**2/4 ) # m^2
    hydr_diameter = clad_diam_out * ( (2*sqrt(3)/pi)*(pin_pitch/clad_diam_out)**2 - 1 ) # m
    avg_velocity = cool_mass_flow/( net_area*density) # m/s

    # adimensional numbers
    Re = density*avg_velocity*hydr_diameter/dyn_viscosity
    Pr = spec_heat*dyn_viscosity/th_cond
    Pe = Pr*Re
    Nu = Nu_cool.subs(Pe_cool,Pe)

    h = Nu*th_cond/hydr_diameter
    output_numbers = array([Re,Pr,Pe,Nu])
    properties = np.array([avg_velocity,net_area,density,dyn_viscosity,spec_heat,th_cond])
    return h, output_numbers,properties


#### ***************** TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS ***************** ####
#
## un giorno modifica pure cladding outer diam (th exp)
def temp_cladding_outer(z,clad_d_out):
    """
    :param z: in m
    :param clad_d_out: in m - Varying parameter! (hot geometry)
    :return: in K
    """
    sec_power_at_z = power_lin_distribution(z) / (pi*clad_d_out)
    temp_coolant_at_z = temp_coolant(z)
    htc,_,_ = heat_transfer_coeff_local(temp_coolant_at_z,clad_d_out) # we need only htc
    temp_clad_out = temp_coolant_at_z + sec_power_at_z/htc
    return temp_clad_out


#### *********** TEMP PROFILE OF CLADDING (INNER) ALONG Z AXIS ************** ####
#
# as we don't know thickness, I'm just using a reasonable value (see design_specifications.py), only later we're going to solve for it...
# anyway to find temperatures changes we only need to set the input and it's done! So quite easy
#
# HP thickness variation due to thermal exp assumed to be constant (to simplify calculations, reasonable approach): eventually demonstrable later anyway....
#
## OSS for hot geometry -> var. of D_clad_outer;
def temp_cladding_inner(z,clad_d_out,clad_thick):
    """
    Give the temperature profile of inner cladding
    :param clad_d_out: Varying parameter! (hot geometry)
    :param z: in [m]
    :param clad_thick: in [m] - UNKNOWN - Initially reasonable value used, then changed...
    :return: in kelvin
    """
    # Initial guess to solve stuff below...
    temp_ci_guess = (600+273.15)

    # writing equation ( stuff ... ... = 0 ) to be solved via minimization by SciPy (see equation_temp_solver method in general_functions.py)
    eqz_1 = temp - temp_cladding_outer(z,clad_d_out) # variable: temp presa da mat_properties...
    eqz_2 = power_lin_distribution(z) * clad_thick / ( pi * (clad_d_outer-2*clad_thick) * clad_thermal_cond ) # clad th cond dep. on temp too!
    res = eqz_1 - eqz_2

    output = equation_temp_solver(res, temp_ci_guess)
    return output


#### ****************** TEMPERATURE PROFILE ALONG OUTER FUEL ******************* ####
#
# NB variazione diam ext del fuel (hot geometry)
# NB considerare pure variazione composizione gas??
#
# NOTE: CONTACT HEAT EXCHANGE NOT IMPLEMENTED, NEITHER RADIATIVE ONE!
##
def temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th):
    """
    DELTA GAP UNKNOWN, it depends on Fuel out and Clad in, which depends itself on Thickness (ASSUMED CONSTANT)
        and Clad out (fixed by project but expanding...)
    :param clad_d_out: Varying parameter! (hot geometry)
    :param z: in meters
    :param fuel_diam_outer: in meters - Varying parameter! (hot geometry)
    :param clad_th: in m - UNKNOWN
    :return: temperature and also delta gap
    """
    temp_fuel_outer_guess = 1000 + 273.15
    delta_gap = (clad_d_out - 2 * clad_th - fuel_diam_outer) / 2
    delta_gap_eff = delta_gap + 10e-6 # m - 10e-6 He as we consider initially 100% He...
    temp_clad_in = temp_cladding_inner(z,clad_d_out,clad_th)

    eqz_1 = temp - temp_clad_in
    eqz_2 = power_lin_distribution(z) * delta_gap_eff / ( pi * fuel_diam_outer * helium_thermal_cond )
    res = eqz_1 - eqz_2
    if delta_gap >= 0:
        out = equation_temp_solver(res, temp_fuel_outer_guess)
    else:
        out = temp_clad_in  # no contact implemented
    return out, delta_gap


#### ******************** TEMPERATURE PROFILE OF OUTER FUEL ALONG Z AXIS ********************** ####
#
# NB variazione diam ext del fuel (hot geometry)
#
## da espandere per bene (void factor, zone restructuring, pu redistri... è una bozza al momento!!)
def temp_fuel_inner(z,clad_d_out,fuel_diam_outer,clad_th):
    """
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: in K
    """
    temp_fuel_inner_guess = 1500 + 273.15

    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th)
    k_fuel = fuel_thermal_cond.subs(x_om,2) # per ora questi valori
    k_fuel = k_fuel.subs(pu_conc,0.2)
    k_fuel = k_fuel.subs(por,0.12)

    ### NB! MODIFICARE QUI PER RISPETTARE EQUAZIONE!!!! SAI COSA
    eqz_1 = temp - temp_fuel_out
    eqz_2 = power_lin_distribution(z) / (4*pi*k_fuel)
    res = eqz_1 - eqz_2
    out = equation_temp_solver(res, temp_fuel_inner_guess)
    return out

#### ******************** TEMPERATURE PROFILE OF OUTER FUEL ALONG RADIUS ********************** ####
#
# HP no azimuthal, angular dependence
# NB variazione diam ext del fuel (hot geometry)
#
## da espandere per bene (void factor, zone restructuring, pu redistri... è una bozza al momento!!)
def temp_fuel_inner_radial(r,z,clad_d_out,fuel_diam_outer,clad_th):
    """
    :param r: in m
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: in K
    """
    temp_fuel_inner_guess = 1500 + 273.15

    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th)
    k_fuel = fuel_thermal_cond.subs(x_om,2) # per ora questi valori
    k_fuel = k_fuel.subs(pu_conc,0.2)
    k_fuel = k_fuel.subs(por,0.12)

    fuel_radius_outer = fuel_diam_outer/2
    eqz_1 = temp - temp_fuel_out
    eqz_2 = ( power_lin_distribution(z) / ( 4*pi*k_fuel ) ) * ( 1 - (r/fuel_radius_outer)**2 )

    res = eqz_1 - eqz_2
    output = equation_temp_solver(res,temp_fuel_inner_guess)
    return output


#### ********************************** HOT GEOMETRY FUNCTION **************************************** ####
def diameter_th_exp_cladding(diam,temperature):
    return diam + diam * alfa_cladding * (temperature - temp_in) ### RIVEDI QUESTA ROBA CHE FORSE MANCA IL DATO!

def diameter_th_exp_fuel(diam,temperature):
    return diam + diam * alfa_fuel * (temperature - temp_in)
