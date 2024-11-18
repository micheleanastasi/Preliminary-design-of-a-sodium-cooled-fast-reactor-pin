"""
#### NOTE ####
c'è ancora da considerare burn up effects
anche fission gases --> K_gas
anche restructuring
"""
from numpy import pi, log
from math import sqrt
from sympy.physics.units import temperature

from general_functions import *
from material_properties import *
from design_specifications import *


#### ************************* TEMP PROFILE OF COOLANT ALONG Z AXIS ************************ ####

def temp_coolant(z):
    """
    Metodo che minimizzi la funzione res affinche si trovi una temperatura che rispetti il bilancio energetico,
    in quanto Cp stesso dipende da temperatura!

    :param z: position choice
    :return: Temp of coolant at point T
    """
    T_0 = 800 + 273.15 # K

    # energy balance: f*int(Cp(T), Tin, T) - int(q', bottom_pin, z) / mass_flow = 0 ... Tout is our goal
    eqz_2 = f_pitch*integral_power_lin_distr(z) / cool_mass_flow  # valore num seconda parte eqz !! NB F PITCH=1/2 !!
    eqz_1 = sy.integrate(cool_spec_heat,temp)
    res = eqz_1 - eqz_1.subs(temp,cool_temp_inlet) - eqz_2 # primo è integtale, secondo intehrale con sostituzione (quindi risolto) e terzo eqz_2

    T_coolant_at_z = equation_solver(res, T_0)
    return T_coolant_at_z



#### ********************* TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS ********************* ####

## CORRELATION FUNCTIONS

# Calculated for each temperature of the coolant along z axis
# HP: pitch btw pins constant, no cross-section of flows!!

def heat_transfer_coefficient(temperature,clad_diam_out):
    """
    Get parameters of coolant such as htc, adim. numbers and avg_velocity: point by point, according to temperature and clad diameter.

    HYPOTHESIS:
        - pitch btw pins constant, no cross-section of flows!!

    :param clad_diam_out: in m - Varying parameter! (hot geometry)
    :param temperature: in [K]
    :return: htc; array with adim. numbers; array with some coolant/geometry properties calculated here
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
    adim_numbers = array([Re, Pr, Pe, Nu])
    properties = np.array([avg_velocity, net_area, density, dyn_viscosity, spec_heat, th_cond])
    return h, adim_numbers, properties


def temp_cladding_outer(z,clad_d_out):
    """
    :param z: in m
    :param clad_d_out: in m - Varying parameter! (hot geometry)
    :return: in K
    """
    sec_power_at_z = power_lin_distribution(z) / (pi*clad_d_out)
    temp_coolant_at_z = temp_coolant(z)
    htc,_,_ = heat_transfer_coefficient(temp_coolant_at_z, clad_d_out)  # we need only htc
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
    eqz_2 = power_lin_distribution(z) * clad_thick / ( pi * (clad_d_out-2*clad_thick) * clad_thermal_cond ) # clad th cond dep. on temp too!
    res = eqz_1 - eqz_2

    output = equation_solver(res, temp_ci_guess)
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
    DELTA GAP UNKNOWN, it depends on Fuel out and Clad in, which depends itself on Thickness (ITS VARIATION ASSUMED NULL)
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
        out = equation_solver(res, temp_fuel_outer_guess)
    else:
        out = temp_clad_in  # no contact implemented yet, so meaningless at the moment
    return out, delta_gap



#### **************************** TEMPERATURE PROFILE OF INNER FUEL  ********************************* ####

## ALONG Z AXIS - SINGLE REGION (NO RESTRUCTURING)
## variabili per conc Pu, porosità ecc
def temp_fuel_inner(z,clad_d_out,fuel_diam_outer,clad_th,x=0,pu=0.29,po=0.12,fv=1):
    """
    Computing inner
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: temp fuel inner in K
    """
    # da espandere (void factor, zone restructuring, pu redistri... è una bozza al momento!!)
    temp_fuel_inner_guess = 1500 + 273.15

    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th) # trattino basso poiché non serve valore delta gap
    k_fuel = fuel_thermal_cond.subs(x_om,x) # per ora questi valori
    k_fuel = k_fuel.subs(pu_conc,pu)
    k_fuel = k_fuel.subs(por,po)

    eqz_1 = temp - temp_fuel_out
    eqz_2 = fv * power_lin_distribution(z) / (4*pi*k_fuel)
    res = eqz_1 - eqz_2
    out = equation_solver(res, temp_fuel_inner_guess)
    return out


## ALONG RADIUS - SINGLE REGION (NO RESTRUCTURING)
def temp_fuel_inner_radial(r,z,clad_d_out,fuel_diam_outer,clad_th,x=0,pu=0.29,po=0.12,r_c=0,temp_clmn=0):
    """
    HP no azimuthal, no angular dependence

    NOTE: from x to po: extraneous parameters to expand functionalities to restructuring (btw they've default, hence no
    restructuring, values)
    :param r: in m
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: in K
    """
    temp_fuel_inner_guess = 1500 + 273.15

    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th)
    k_fuel = fuel_thermal_cond.subs(x_om,x) # per ora questi valori
    k_fuel = k_fuel.subs(pu_conc,pu)
    k_fuel = k_fuel.subs(por,po)

    fuel_radius_outer = fuel_diam_outer/2
    if temp_clmn == 0:
        eqz_1 = temp - temp_fuel_out
    else:
        eqz_1 = temp - temp_clmn
    eqz_2 = ( power_lin_distribution(z) / ( 4*pi*k_fuel ) ) * ( 1 - ((r)/fuel_radius_outer-r_c)**2 ) # sensible?? r - rclmn

    res = eqz_1 - eqz_2
    output = equation_solver(res, temp_fuel_inner_guess)
    return output



#### ********************************** DIAMETER - HOT GEOMETRY FUNCTION ***************************************** ####

def diameter_th_exp_cladding(diam,temperature):
    clad_exp = clad_eps_th.subs(temp,temperature)
   # return diam + clad_exp*diam
    return diam + diam * alfa_clad * (temperature - temp_in)

def diameter_th_exp_fuel(diam,temperature):
    return diam + diam * alfa_fuel * (temperature - temp_in)



#### *********************************** ITERATION - HOT GEOMETRY FUNCTION **************************************** ####
def hot_geometry_iteration(z, clad_d_out_0, fuel_d_out_0, clad_thick_0,print_status=True):
    """
    Iterative calculation in order to get temperatures, gap, final fuel out and clad out diameters, other properties
    (such as htc, net area, adim. numbers and so on...) in a hot geometry model
    NOTE! FOR RESTRUCTURING USEFUL UNTIL FUEL OUTER!
    :param z: in m
    :param clad_d_out_0: cold geo diameter, in m
    :param fuel_d_out_0: cold geo diameter, in m
    :param clad_thick_0: supposed cladding thickness to be used, in m
    :return: cold temp, hot temp, minimun gap along the pin, generic properties @ hot geo, final clad outer diam, final fuel outer diam
    """
    tol = 10e-3
    temp_array = np.zeros(5) # 0: temp coolant - 1: temp clad out - 2: temp clad in - 3: temp fuel out - 4: temp fuel in

    # cold geometry computing - to be used as first iteration data
    temp_array[0] = temp_coolant(z)
    temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
    temp_array[2] = temp_cladding_inner(z, clad_d_out_0, clad_thick_0)
    temp_array[3],_ = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thick_0)
    temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thick_0)
    old_temp = temp_array.copy() # to give as output the cold geo temps too

    while True:
        prec_temp_array = temp_array.copy() # to be used to evaluate when exiting from the while below...

        # considerare sempre espansione rispetto al diametro INIZIALE
        clad_d_out_0 = diameter_th_exp_cladding(clad_d_outer, prec_temp_array[1])  # with temp clad outer
        fuel_d_out_0 = diameter_th_exp_fuel(fuel_d_outer, prec_temp_array[3])  # with temp fuel outer

        # hot geo computing
        temp_array[0] = temp_coolant(z)
        temp_array[1] = temp_cladding_outer(z, clad_d_out_0)
        yy_htc_loc, yy_adim_num_cool, yy_cool_loc_prop = heat_transfer_coefficient(temp_array[0], clad_d_out_0)  # new
        temp_array[2] = temp_cladding_inner(z, clad_d_out_0, clad_thick_0)
        temp_array[3], delta_gap = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thick_0)
        temp_array[4] = temp_fuel_inner(z, clad_d_out_0, fuel_d_out_0, clad_thick_0)

        if np.abs(prec_temp_array[4] - temp_array[4]) < tol:  # va bene così (?)
            break
    other = np.array(list([yy_htc_loc]) + list(yy_adim_num_cool) + list(yy_cool_loc_prop))

    if print_status: # optional "progress bar" print (see input boolean)
        print(f"Hot geo completed at {np.round(100 * z / 0.85, 2)}% (Position: {np.round(z, 2)} m) - Temp fuel "
              f"inner: HOT:{np.round(temp_array[4], 2)} K, COLD:{np.round(old_temp[4], 2)} K - Min. gap:{np.round(delta_gap * 1e6,3)} um")

    return old_temp, temp_array, delta_gap, other, clad_d_out_0, fuel_d_out_0



#### ********************************************** RESTRUCTURING ************************************************* ####

## preparatory functions
def void_factor(rad_fv,rad_fo):
    x = (rad_fo/rad_fv)**2
    return 1 - ( log(x) )/(x - 1)


def get_R_from_temp(z,r_fuel_out,temperature,T_fuel_in,T_fuel_out):
    """
    Useful to get radius corresponding to a certain temp in the mono-region pellet model
    Eventually useful to divide in regions
    :param z:
    :param r_void:
    :param r_fuel_out:
    :param temperature:
    :param T_fuel_out:
    :return:
    """
    check = T_fuel_in > temperature
    if check:
        k_fuel = fuel_thermal_cond.subs(x_om,2) # mettere variabili
        k_fuel = k_fuel.subs(pu_conc,0.29)
        k_fuel = k_fuel.subs(por,0.12)
        k_fuel = k_fuel.subs(temp,temperature)
        output = r_fuel_out*sqrt(1 - (temperature-T_fuel_out) * (4 * pi * k_fuel) / ( power_lin_distribution(z) ))
    else:
        output = 0
    return output


def radius_void_get(r_clmn,porosity):
    """
    to get R void
    :param r_clmn:
    :param porosity:
    :return:
    """
    output = sqrt( porosity ) * r_clmn # get R void
    return output


def fuel_region_temp_radial(z,r):

    return None