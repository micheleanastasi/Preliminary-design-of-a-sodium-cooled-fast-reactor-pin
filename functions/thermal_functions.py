"""
Collection of functions aimed at computing thermal properties
"""

import numpy as np
from mpmath import arange
from numpy import pi,round
from scipy.integrate import quad as integral
from math import sqrt,log
import os

from sympy.vector import Divergence

from functions.general_functions import *
from functions.general_properties import *

#### ************************************************************************************************************** ####
#### ***************************************** PREPARATORY FUNCTIONS ********************************************** ####
#### ************************************************************************************************************** ####

## TEMP PROFILE OF COOLANT ALONG Z AXIS
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

    T_coolant_at_z = sy_equation_solver(res, T_0)
    return T_coolant_at_z



## CORRELATION FUNCTIONS FOR TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS
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


## TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS
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



## TEMP PROFILE OF CLADDING (INNER) ALONG Z AXIS
def temp_cladding_inner(z,clad_d_out,clad_thick):
    """
    Give the temperature profile of inner cladding

    NOTE:
    - as we don't know thickness,  just using a reasonable value (see general_properties.py), only later we're going
      to solve for it... anyway to find temperatures changes we only need to set the input and it's done! So quite easy

    HP:
    - thickness variation due to thermal exp assumed to be constant (to simplify calculations, reasonable approach), but
      actually its value may change by implementing swelling phenomenum (of course see successive functions...)

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

    output = sy_equation_solver(res, temp_ci_guess)
    return output



## TEMPERATURE PROFILE ALONG OUTER FUEL - DA COMPLETARE
def temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th,burnup,contact_pressure):
    """
    NOTE:
    - CONTACT HEAT EXCHANGE NOT YET IMPLEMENTED, NEITHER RADIATIVE ONE!
    - DELTA GAP UNKNOWN firstly, since it depends on Fuel out and Clad in, which depends itself on Thickness (ITS
        VARIATION ASSUMED NULL) and Clad out (fixed by project but expanding...)

    :param clad_d_out: Varying parameter! (hot geometry)
    :param z: in meters
    :param fuel_diam_outer: in meters - Varying parameter! (hot geometry)
    :param clad_th: in m - UNKNOWN
    :return: temperature and also delta gap
    """
    temp_fuel_outer_guess = 1000 + 273.15

    delta_gap = (clad_d_out - 2 * clad_th - fuel_diam_outer) / 2
    delta_gap_eff = delta_gap + 10e-6 # m - 10e-6 He as we consider initially 100% He... RVD CIO!!!! cambio % he
    roughness=10e-6 #[m]

    temp_clad_in = temp_cladding_inner(z,clad_d_out,clad_th)
    m_xe, m_kr, m_he = fg_prod(burnup)
    gap_k = k_th_gas(temp_clad_in,m_he,m_xe,m_kr) # (HPCONS)

    eqz_1 = temp - temp_clad_in
    if delta_gap > 0:
        eqz_2 = power_lin_distribution(z) * delta_gap_eff / (pi * fuel_diam_outer * gap_k)
        res = eqz_1 - eqz_2
        out = sy_equation_solver(res, temp_fuel_outer_guess)
    else:
        #out = temp_clad_in + 400 # approx for contact
        # worst case: very few contact points (edges!) then for the rest gap is like roughness size (order of magn: micron)
        #eqz_2 = power_lin_distribution(z) * 5e-6 / (pi * fuel_diam_outer * gap_k)
        k_f=k_th_fuel(temp_fuel_outer_guess,burnup)
        k_c=clad_thermal_cond.subs(temp,temp_clad_in)
        gap_h=gap_k/roughness+fg_C*2*k_f*k_c/(k_f+k_c)*contact_pressure/fuel_hardness/roughness**0.5
        eqz_2 = power_lin_distribution(z) / (pi * fuel_diam_outer * gap_h)

        res = eqz_1 - eqz_2
        out = sy_equation_solver(res, temp_fuel_outer_guess)
    return float(out), delta_gap, gap_k



## TEMPERATURE PROFILE OF INNER FUEL - SINGLE REGION (NO RESTRUCTURING, see later)
def temp_fuel_max(z,clad_d_out,fuel_diam_outer,clad_th,burnup,contact_pressure,fv=1):
    """
    Computing inner according to z position
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: temp fuel inner in K
    """
    temp_0 = 1500 + 273.15

    temp_fuel_out,_,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th,burnup,contact_pressure) # trattino basso poiché non serve valore delta gap
    res = lambda t : t - temp_fuel_out - ((t - temp_fuel_out) * fv * power_lin_distribution(z) / (4 * pi * integral(k_th_fuel, temp_fuel_out, t, args=(burnup,))[0]))
   # out  = temp_fuel_out + (fv * power_lin_distribution(z) / (4 * pi * k_th_fuel(temp_fuel_out,burnup)) )

    out = fun_equation_solver(res, temp_0)
    return out

## PLOT RADIAL
def temp_fuel_radial(rr,fuel_temp_out,fuel_temp_max,fuel_r_out,fuel_r_void):

    a=(fuel_temp_out-fuel_temp_max)*(fuel_r_out-fuel_r_void)**-2
    b=-2*fuel_r_void*a
    c=fuel_temp_max+a*fuel_r_void**2
    tempfuelradial= a * rr ** 2 + b * rr + c

    tempfuelradial[rr < fuel_r_void]=max(tempfuelradial)

    return tempfuelradial



#### ******************************************** GEOMETRIES EXPANSIONS ******************************************* ####

def diameter_th_exp_cladding(diam,t_max):
    """
    NOTE:
    - should be assumed as temp the one related to clad inner, to be conservative!
    """
    clad_exp = clad_eps_th.subs(temp,float(t_max))
    return diam + clad_exp*diam


def diameter_th_exp_fuel(z,diam,t_max, t_min, burnup,rad_void):
    """
    CONSIDERING THERMAL EXPANSION
    By equation: new_rad = old_rad + alfa * Int(A-Br^2, btw diam_Void/zero and old rad)
    OSS burn up in GWd/ton(HM)
    """
    radius = diam/2
    k_integ = integral(k_th_fuel, t_min, t_max, args=(burnup,))[0]
    A = t_min + (t_max-t_min)*power_lin_distribution(z)/( 4*pi*k_integ )
    B = (t_max-t_min)*power_lin_distribution(z)/( 4*pi*k_integ*radius**2 )

    fun = lambda r : A - B * r**2
    integ = integral(fun,rad_void,diam/2)[0]
    out = radius + alfa_fuel * integ

    #out = ( radius + alfa_fuel * ( A*radius - (B/3)*radius**3 ) )

    return out * 2


def length_th_exp_cladding(leng,t_max):
    """
    To be used to take into account increase of volume for pressure
    """
    return leng + leng*clad_eps_th.subs(temp,t_max)

def length_th_exp_fuel(leng,t_max):
    """
    CONSIDERING THERMAL
    OSS burn up in GWd/ton(HM)
    """
    out = ( leng + leng*alfa_fuel* ( t_max - temp_in) ) # (HP CONS)
    return out * 2




#### ********************************************** RESTRUCTURING ************************************************* ####

def void_factor(rad_fv,rad_fo):
    x = (rad_fo/rad_fv)**2
    return 1 - ( log(x) )/(x - 1)

def get_R_from_temp(z,d_fuel_out,temperature,T_fuel_in,T_fuel_out,burnup):
    """
    Useful to get radius corresponding to a certain temp in the mono-region pellet model
    Eventually useful to divide in regions
    """
    r_fuel_out = d_fuel_out/2
    check = T_fuel_in > temperature # does it happen?
    if check:
        k_fuel = integral(k_th_fuel, T_fuel_out, temperature, args=(burnup,))[0] / (temperature - T_fuel_out)
        output = r_fuel_out*sqrt(1 - (temperature-T_fuel_out) * (4 * pi * k_fuel) / ( power_lin_distribution(z) ))
    else:
        output = 0
    return output


def radius_void_get(r_clmn,porosity):
    output = sqrt( porosity ) * r_clmn # get R void
    return output


def fuel_restructuring(z,temp_fuel_out,temp_fuel_in,diam_fuel_out,diam_clad_out,burnup,print_stuff=False,contact_pressure=0):
    """
    HP CONS:
    - using fv but neglecting beneficial increase of k_Fuel in col. region...
    - no equi-axial region designed

    OSS:
    - no changing radius with time (burnup), hence generating for low burnup (about 1), then saving them to be used in b-up higher levels, otherwise we would get lower diameters
        even if the initial radius should not change so dramatically and this would be unphysical, despite the fact that it would be too conservative
    - reduction of void volume (here normalized per height, so area!) to make up for the contact with the cladding (HP SEMPL): pellet, ductile inside (NB NDTT...) and cracked, can
        be thought as closing the void region to accomodate swelling expansion

    GET:
    - new diameter of fuel by only ONE iteration (like hot geo iteration)
    - old diamtere of fuel
    - radius of columnar region
    - radius of void region
    - temperature @ void, hence new maximum: got using temp_fuel_max + void factor
    - old array of temperatures
    - interference if null void region eventuall
    """
    # creating radius_void to be used for other (higher) level of burn up, as said above
    if 0.8 <= burnup < 1.3:
        radius_clmn = get_R_from_temp(z, diam_fuel_out, fuel_temp_clmn, temp_fuel_in, temp_fuel_out, burnup)
        radius_void = radius_void_get(radius_clmn, poro_asf)
        np.save(os.path.join("restr_saves",f"radius_clmn_{z}_burnup.npy"),radius_clmn)
        np.save(os.path.join("restr_saves",f"radius_void_{z}_burnup.npy"),radius_void)
    elif burnup >= 1.5:
        try:
            radius_clmn = np.load(os.path.join("restr_saves",f"radius_clmn_{z}_burnup.npy"))
            radius_void = np.load(os.path.join("restr_saves",f"radius_void_{z}_burnup.npy"))

        except FileNotFoundError:
            print(f"ERROR! FILE REGARDING LOW BURNUP VOID AND CLMN RADIUS NOT FOUND! First re-run with 0.8 <= burnup < 1.3, following data and"
              f"results can be inaccurate, since null value is assigned to the two radius")
            radius_clmn = 0
            radius_void = 0
    else:
        radius_clmn = 0
        radius_void = 0

    old = temp_fuel_max(z, diam_clad_out, diam_fuel_out, clad_thickness_0,burnup,contact_pressure)
    old_diam = diam_fuel_out

    if radius_clmn != 0 and radius_void != 0:
        voidFactor = void_factor(radius_void,diam_fuel_out/2)
       # temp_void = temp_fuel_max(z,diam_clad_out,diam_fuel_out,clad_thickness_0,burnup,voidFactor)

        temp_void = temp_fuel_out + (temp_fuel_in - temp_fuel_out) * voidFactor

        ## hot geo - only one "iteration"
        #diam_fuel_out = diameter_th_exp_fuel(z, fuel_d_outer, temp_void, temp_fuel_out,burnup,)

        if print_stuff:
            print(f"***** ITERATION *****")
            print(f"\n>DATA ABOUT RESTRUCTURING:")
            print(f"Temp void: {round(float(temp_void), 2)} K - Temp old: {np.round(temp_fuel_in, 2)} K")
            print(f"Radius of clmn: {round(float(radius_clmn),6)} m, radius of void: {round(float(radius_void),6)} m")
            print(f"*********************")

    else:
        temp_void = temp_fuel_in
        if print_stuff:
            print("\n>NO RESTRUCTURING HAPPENED")

    return diam_fuel_out, old_diam, radius_clmn, radius_void, float(temp_void), old




#### **************************************************** GAS IN GAP ********************************************** ####

def gap_vol_cold():
    """
    Used to calculate initial moles of He in gap (as you'll see later, we'll also take into account the plenm extra vol)
    """
    vol_1 = pin_column_height * 0.25*pi*clad_d_outer**2
    vol_2 = pin_column_height * 0.25*pi*fuel_d_outer**2 # HP coincidenza altezze
    return vol_1 - vol_2


def gap_vol_hot(fuel_d_out_array,clad_d_out_array,gap):
    """
    NB Arrays as input!
    """
    num = len(fuel_d_out_array)
    unit = pin_column_height/num    # unit of length (discretization)
    vol = 0     # init volume
    clad_diameter_in = clad_d_out_array - 2*clad_thickness_0

    if np.sum(gap) <= 0:    #obviously no gap volume...
        return 0
    else:
        for i in range(0,num):
            # condition to avoid summing up a negative volume
            if clad_diameter_in[i] > fuel_d_out_array[i]:
                vol = vol + unit * 0.25 * pi * (clad_diameter_in[i] ** 2 - fuel_d_out_array[i] ** 2)

    return vol

def fg_prod(burnup):
    """
    Burn up in GWd/ton!!

    HP: FGR = 100 % (conservative/worst case and also easier)
    """
    FGR=1 # fraction of produced gases that are released

    N_av = 6.022e23 # atoms/mol
    y_xe = 0.27
    y_kr = 0.03
    y_he = 0.022

    E_fiss = 200e6 * 1.6e-19 # J
    # q_third = power_lin_max/ (pi*0.25*fuel_d_outer**2) # W/m^3

    #   consider avg power for gas production
    # print(f'avg power {power_lin_avg}')
    q_third = power_lin_avg/ (pi*0.25*fuel_d_outer**2) # W/m^3
    vol = pin_column_height * (pi*0.25*fuel_d_outer**2) # m^3
    F_prime = q_third / E_fiss # fiss/s/m^3

    mass = fuel_density * vol/1000 *(239.77/271.77)  # ton (heavy metal), last factor is MM_metal/MM_MOX
    q = power_lin_max * pin_column_height / 1e9 # GW
    time = burnup * 86400 * mass/ q # seconds

    # print(f'time : {time}')

    mol_xe = FGR*( y_xe * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    mol_kr = FGR*( y_kr * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    mol_he = FGR*( y_he * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    return mol_xe*time, mol_kr*time, mol_he*time



#### ************************************************************************************************************** ####
#### ******************************************* PRINCIPAL FUNCTIONS ********************************************** ####
#### ************************************************************************************************************** ####

def thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0,bup,contact_pressure):
    """
    Use it to get principal points temperatures: if used once, then you'll get cold geo temps, hence you should embed it
    in an iterative cycle properly set up in order to implement and get hot geometry temps.
    Also giving delta gap measures
    """
    output = np.zeros(5) # 0: temp coolant - 1: temp clad out - 2: temp clad in - 3: temp fuel out - 4: temp fuel in

    output[0] = temp_coolant(z)
    output[1] = temp_cladding_outer(z, clad_d_out_0)
    output[2] = temp_cladding_inner(z, clad_d_out_0, clad_thick_0)
    output[3], d_gap,_ = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thick_0, bup,contact_pressure)
    output[4] = temp_fuel_max(z, clad_d_out_0, fuel_d_out_0, clad_thick_0, bup,contact_pressure)
    return output, d_gap


def hot_geometry_general(z, clad_d_out_0, fuel_d_out_0, clad_thick_0,bup,print_status=True,RestructOn=True):
    """
    This is one of the more important functions implemented here:
    Iterative calculation in order to get temperatures, gap, final fuel out and clad out diameters, other properties
    (such as htc, net area, adim. numbers and so on...) in a hot geometry model


    NOTE:
    - Also accounting for fuel and cladding swelling
    - FOR RESTRUCTURING USEFUL UNTIL FUEL OUTER! (then substitution of max fuel temp by restr function output)
    - NEW! Implementation of burnup (our "time") eventually available, although by default set at zero (see above)
    - set RestructOn to take into account restructuring effects (however not computed if under a centain threshold)

    OSS
    - with higher burnup we would have lower levels of restructuring (as noted a posteriori, temperatures generally decrease): then I do not consider
        the radius values (clmn, void) for burnup = 0 but the current bup value of the function, in order to minimize the beneficial effecys (HP CONS) and
        especially to simplify the code!!
    - about order of properties:
      'HTC [W/m^2/K]','Re','Pr','Pe','Nu', 'avg_velocity [m/s]','net_area [m^2]','density [kg/m^3]','dyn_viscosity [Pa*s]','spec_heat [J/kg/K]', 'th_cond [W/m/K]']

    :param z: in m
    :param clad_d_out_0: cold geo diameter, in m
    :param fuel_d_out_0: cold geo diameter, in m
    :param clad_thick_0: supposed cladding thickness to be used, in m
    :param bup: Burnup in GWd/ton
    :return: cold temp, hot temp, minimun gap along the pin, generic properties @ hot geo, final clad outer diam, final fuel outer diam
    """

    # VARIABLE INITIALIZATION
    isOkayColdGeo = False # used to consider only cold geo if delta gap < 0 (HP CONS)
    isFirstTime = True
    isSurelyNegative = False
    radClmn = 0
    interf = 0
    conv = 0
    old_clad = clad_d_out_0
    old_fuel = fuel_d_out_0

    # SETTING
    tol = 5 # kelvin - tolerance set for following iteration phase (see later in this def)

    # to simplify the code, starting from a guess value of rad void (then correction by iteration)
    if bup > 1.5:
        radVoidGuess = 0.5e-3  # m
    else:
        radVoidGuess = 0

    # 0: temp coolant - 1: temp clad out - 2: temp clad in - 3: temp fuel out - 4: temp fuel in
    # first data useful for swelling, cold geo computing
    temp_array,old_gap = thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0, bup,contact_pressure=0)
    old_temp = temp_array.copy() # to give as output the cold geo temps too


    # variation of external diameters due to swelling
    # for cladding increase of thickness! (not more than about 4 %)
    fuel_d_outer = swelling_fuel(z,bup, fuel_d_out_0)
    clad_d_outer = swelling_clad(z, bup, temp_array[2], clad_d_out_0)
    delta_gap = old_gap


    ## CONSOLE PRINT SET UP here
    if print_status:
        print(f">STARTING COMPUTING @ {np.round(100 * z / 0.85, 2)}% (POSITION: {np.round(z, 2)} m)")


    # START OF ITERATIVE COMPUTING OF TEMPERATURES #
    iter = 0 # iteration counting
    while True:
        iter += 1
        prec_temp_array = temp_array.copy()  # to be used to evaluate when exiting from the while below...
        prec_conv = float(conv)
        prec_gap = delta_gap

        # loading r_void already computed if burnup higher than the below threshold (r_void used to compute contact pressure)
        if bup < 1.3:
            fuel_r_void = np.zeros(z.shape)
        else:
            fuel_r_void = np.load(os.path.join("restr_saves",f"radius_void_{z}_burnup.npy"))

        # always consider expansion w.r.t. initial geometry (clad_d_outer, fuel_d_outer variables)
        clad_d_out_0 = diameter_th_exp_cladding(clad_d_outer, (prec_temp_array[2]))
        fuel_d_out_0 = diameter_th_exp_fuel(z, fuel_d_outer, prec_temp_array[4], prec_temp_array[3], bup, radVoidGuess)
        clad_d_in_0 = clad_d_out_0 - 2 * clad_thick_0

        # computing contact pressure
        fuel_E = fuel_Young_modulus(prec_temp_array[3],poro_asf)
        clad_E = clad_Young_modulus(prec_temp_array[2])
        clad_nu = clad_Poisson_ratio(prec_temp_array[2])
        fc_contact_pressure = contact_pressure(fuel_d_out_0 / 2, clad_d_in_0 / 2, clad_d_out_0 / 2, fuel_nu, fuel_E, clad_nu, clad_E, fuel_r_void)

        # hot geo computing, starting by new geometries above
        temp_array, delta_gap = thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0, bup,fc_contact_pressure)

        # considering now restructuring
        rst_threshold = 0.1  # GWd/ton - after about 10 days...
        if RestructOn and bup >= rst_threshold:  # put here since not so variation of fuel diameter following restructuring... (HP SEMPL)
            # giving diameter after a single cycle of iter following restructuring and new temperature (void factor)
            # moreover overwriting void diameter guess value (if first iteration)
            _, _, radClmn, radVoidGuess, temp_array[4], oldMaxTemp = fuel_restructuring(z, temp_array[3], temp_array[4], fuel_d_out_0, clad_d_out_0, bup, print_stuff=False)


        # ITERATION SAFETY
        conv = np.abs(prec_temp_array[4] - temp_array[4])
        #checking to avoid code divergence (or blocked convergence) happening when gap size is negative, running only another time though
        #happening if gap value, got from iteration, oscillating around zero (from pos to neg and vice versa) for every iter
        divergenceSafety = np.abs(conv) > np.abs(prec_conv) and isFirstTime==False
        notConvergenceSafety = (delta_gap < 0 and iter > 40) or delta_gap < 0 and np.abs(prec_conv-conv)<1
        gapAlmostZero = np.abs(delta_gap) < 0.5e-6 and np.abs(prec_conv-conv)<1

        if divergenceSafety: # all sensitive conditions to avoid divergence...
            isOkayColdGeo = True
            fuel_d_out_0 = clad_d_out_0 - 2 * clad_thick_0 - 2 * delta_gap # used to reset to perform further calcs... NB clad temps not influenced by what's inside,so same value... then delta gap null!
            fuel_d_out_0 = float(fuel_d_out_0) # code would crash without this conversion :(
            fc_contact_pressure = contact_pressure(fuel_d_out_0 / 2, clad_d_in_0 / 2, clad_d_out_0 / 2, fuel_nu, fuel_E,clad_nu, clad_E, fuel_r_void)
            temp_array, _ = thermal_computing(z, clad_d_out_0, fuel_d_out_0, clad_thick_0, bup,fc_contact_pressure) # de facto now gap equal to zero for thermal calc but not for mech

        if notConvergenceSafety: # if oscillating around zero - approx: linearizing around zero
            gapAlmostZero = True
            for i in range(0,len(temp_array)):
                temp_array[i] = (temp_array[i]+prec_temp_array[i])/2
            delta_gap = (delta_gap+prec_gap)/2
            fuel_d_out_0 = clad_d_out_0 - 2 * clad_thick_0 - 2 * delta_gap  # used to reset to perform further calcs... NB clad temps not influenced by what's inside,so same value... then delta gap null!
            fuel_d_out_0 = float(fuel_d_out_0)  # code would crash without this conversion :(
            fc_contact_pressure = contact_pressure(fuel_d_out_0 / 2, clad_d_in_0 / 2, clad_d_out_0 / 2, fuel_nu, fuel_E,clad_nu, clad_E, fuel_r_void)

        if np.abs(prec_conv-conv)<1 and delta_gap>0:    # to speed up code if too low to converge even with open gap...
            gapAlmostZero = True
        if iter == 20:
            print(f"Do not worry, code still running")
        if iter > 80:
            print(f"ERROR! DIVERGENCE!")
            break

        # DEBUG INFO
        if True:
            print(f"**convergence: {prec_temp_array[4] - temp_array[4]}")
            print(f"*temp inner {temp_array[4]}")
            print(f"*gap: {delta_gap}")


        # TOLLERANCE CHECK
        if  conv < tol or isOkayColdGeo or gapAlmostZero:
            break
        isFirstTime = False
    # END OF ITERATIVE COMPUTING OF TEMPERATURES #

    # if gap negative assessed, then get as output fuel diam out == fuel clad in, whereas gap negative used for thermal computing (see fuel_outer_temp function)
    if delta_gap < 0 and iter <= 20:
        interf = -delta_gap
        fuel_d_out_0 = clad_d_out_0 - 2 * clad_thick_0  # contact! NB clad temps not influenced by what's inside,so same value... then delta gap null!
        fuel_d_out_0 = float(fuel_d_out_0)  # code would crash without this conversion :(


    # preparing other stuff to be exported by the function
    yy_htc_loc, yy_adim_num_cool, yy_cool_loc_prop = heat_transfer_coefficient(temp_array[0], clad_d_out_0)  # new
    other = np.array(list([yy_htc_loc]) + list(yy_adim_num_cool) + list(yy_cool_loc_prop))



    ## CONSOLE PRINT SET UP here
    if print_status: # optional "progress bar" print (see input boolean)
        print(f"*************** RESULTS @ {np.round(100 * z / 0.85, 2)}% (POSITION: {np.round(z, 2)} m) *************************************************")
        print(f"\n>GEO VARIATION DUE TO SWELLING:\nFuel diameter outer @ cold: {round(float(old_fuel),8)} m --> Now: {round(float(fuel_d_outer),8)} m")
        print(f"Cladding d outer @ cold: {old_clad} m --> Now: {round(float(clad_d_outer),8)} m")
        if radClmn != 0: # meaning that restr happened...
            print(f"\n>RESTRUCTURING HAPPENED AT LOW BURNUP (current vs temp if void region not present)")
            print(f"Temp void: {round(float(temp_array[4]), 2)} K - Temp old: {np.round(oldMaxTemp, 2)} K")
            print(f"Radius of clmn: {round(float(radClmn), 6)} m, radius of void: {round(float(radVoidGuess), 6)} m")
        print(f"\n>FINAL DATA:")
        if divergenceSafety:
            print(f"WARNING! Divergence of code, hence potential positive thermal feedback!")
        print(f"Fuel temp inner: HOT:{np.round(temp_array[4], 2)} K\nOld gap: {round(float(old_gap*1e6),4)} um --> New gap: {round(float(delta_gap*1e6),2)} um"
              f" ({round(float(100 * delta_gap / initial_delta_gap),2)}%)")
        if delta_gap < 0:
            print(f"Interference: {round(float(-interf)*1e6,2)} um")
        print(f"Fuel diam ext: {round(float(fuel_d_out_0*1e3),5)} mm")
        print(f"Clad diam ext: {round(float(clad_d_out_0*1e3),5)} m")
        print(f"\n********************************************************************************************************************************************************\n\n\n\n")

    return old_temp, temp_array, delta_gap, other, clad_d_out_0, fuel_d_out_0, fc_contact_pressure



def pressure_gap_calc(volume_gap,temp_clad_array,burnup,plenum_vol=0,plenum_clad_d_in=0,temp_plenum=0,print_stuff=False):
    """
    Partial pressure computing and then summing them up
    HPs:
    - total release of gas (and also considering swelling, just see other section of the code)
    - temp plenum == coolant one
    - equation of perfect gas
    - not considering axial thermal expansion

    OSS:
    - temp_clad_array is matrix of temp along r and z respectively (rows and columns)
    - plenum_clad_d_in @ temp_env, then expansion will be computed here!
    """
    R = 8.31446 # J/K/mol

    # computing mean temp of gap (hence gas) firstly along radius and then along axis -- only fuel rod section
    mean_temp = np.mean(temp_clad_array,axis=1)
    mean_temp = np.mean(mean_temp)

    # computing initial moles of gas according to HW data
    mol_he_in = (fill_gas_press_in * (plenum_vol + gap_vol_cold())) / (R * fill_gas_temp_in)
    # computing current amount of moles
    mol_xe = fg_prod(burnup)[0]
    mol_kr = fg_prod(burnup)[1]
    mol_he = fg_prod(burnup)[2] + mol_he_in

    # computing a mean temp weighted according to the volumes of fuel pin section and plenum
    vol = volume_gap + plenum_vol
    mean_temp = (volume_gap * mean_temp + plenum_vol * temp_plenum) / (volume_gap + plenum_vol)

    plenum_length = 0   # variable initialization
    plenum_clad_d_in = diameter_th_exp_cladding(plenum_clad_d_in,temp_plenum) # considering also th expansion...
    if plenum_clad_d_in != 0: # condition to simply avoid zero at denominator :D
        plenum_length = plenum_vol / (0.25 * pi * plenum_clad_d_in ** 2) # m

    press_xe = (mol_xe * R * mean_temp) / vol
    press_kr = (mol_kr * R * mean_temp) / vol
    press_he = (mol_he * R * mean_temp) / vol
    total_pressure = press_he + press_xe + press_kr

    if print_stuff:
        print(f"\n************** DATA ABOUT PRESSURE @ {burnup} GWd/ton **************")
        print(f"Pressure: {round(float(total_pressure*1e-6),2)} MPa -- Plenum length: {round(float(plenum_length),2)} m\n")
        print(f"Volume tot: {round(float(vol),9)} m^3 -- Volume gap: {round(float(volume_gap),9)} m^3 -- Volume plenum: {round(float(plenum_vol),9)} m^3")

    return total_pressure, plenum_length  # Pa - m