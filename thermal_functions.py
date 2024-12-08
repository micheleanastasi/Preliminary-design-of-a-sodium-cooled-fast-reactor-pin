"""
Collection of functions aimed at computing thermal properties
"""
from mpmath import arange
from numpy import pi,round
from scipy.integrate import quad as integral
from math import sqrt,log

from scipy.special import y0_zeros

from general_functions import *
from general_properties import *

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



## TEMPERATURE PROFILE ALONG OUTER FUEL
def temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th,burnup=0):
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

    temp_clad_in = temp_cladding_inner(z,clad_d_out,clad_th)
    #gap_k = helium_thermal_cond.subs(temp,temp_clad_in)
    m_xe, m_kr, m_he = fg_prod(burnup=burnup)
   # gap_k = th_gap(temp_clad_in,x_he=m_he,x_kr=m_kr,x_xe=m_xe) # (HPCONS) - AGGIORNAREEEEEEEEEEEEEEEEEEE
    gap_k = k_th_gas(temp_clad_in)

    eqz_1 = temp - temp_clad_in
    #eqz_2 = power_lin_distribution(z) * delta_gap_eff / ( pi * fuel_diam_outer * helium_thermal_cond )
    eqz_2 = power_lin_distribution(z) * delta_gap_eff / (pi * fuel_diam_outer * gap_k)
    res = eqz_1 - eqz_2
    if delta_gap >= 0:
        out = sy_equation_solver(res, temp_fuel_outer_guess)
    else:
        out = temp_clad_in + 100 # no contact implemented yet, so meaningless at the moment
    return out, delta_gap




## TEMPERATURE PROFILE OF INNER FUEL - SINGLE REGION (NO RESTRUCTURING, see later)
def temp_fuel_max(z,clad_d_out,fuel_diam_outer,clad_th,fv=1,burnup=0):
    """
    Computing inner according to z position
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: temp fuel inner in K
    """
    temp_0 = 1500 + 273.15

    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th) # trattino basso poiché non serve valore delta gap
    res = lambda t : t - temp_fuel_out - ((t - temp_fuel_out) * fv * power_lin_distribution(z) / (4 * pi * integral(k_th_fuel, temp_fuel_out, t, args=(burnup,))[0]))

    out = fun_equation_solver(res, temp_0)
    return out


## ALONG RADIUS - SINGLE REGION (NO RESTRUCTURING, see later)
def temp_fuel_inner_radial(r,z,clad_d_out,fuel_diam_outer,clad_th,burnup=0):
    """
    HP:
    - no azimuthal
    - no angular dependence
    - q''' uniform

    :param r: in m
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN --- all these parameters used to calculate delta gap
    :return: in K
    """
    temp_0 = 1500 + 273.15
    fuel_radius_outer = fuel_diam_outer/2
    temp_fuel_out,_ = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th) # trattino basso poiché non serve valore delta gap
    res = lambda t : t - temp_fuel_out - ((t - temp_fuel_out) * power_lin_distribution(z) / (4 * pi * integral(k_th_fuel, temp_fuel_out, t, args=(burnup,))[0])) * (1 - (r / fuel_radius_outer) ** 2)

    output = fun_equation_solver(res, temp_0)
    return output



#### ******************************************** GEOMETRIES EXPANSIONS ******************************************* ####

def diameter_th_exp_cladding(diam,t_max):
    """
    NOTE:
    - should be assumed as temp the one related to clad inner, to be conservative!
    """
    clad_exp = clad_eps_th.subs(temp,t_max)
    return diam + clad_exp*diam


def diameter_th_exp_fuel(z,diam,t_max, t_min, burnup=0):
    """
    CONSIDERING THERMAL EXPANSION
    By equation: new_rad = old_rad + alfa * Int(A-Br^2, btw 0 and old rad)
    OSS burn up in GWd/ton(HM)
    """
    radius = diam/2
    k_integ = integral(k_th_fuel, t_min, t_max, args=(burnup,))[0]
    A = t_min + (t_max-t_min)*power_lin_distribution(z)/( 4*pi*k_integ )
    B = (t_max-t_min)*power_lin_distribution(z)/( 4*pi*k_integ*radius**2 )
    out = ( radius + alfa_fuel * ( A*radius - (B/3)*radius**3 ) )
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

def get_R_from_temp(z,d_fuel_out,temperature,T_fuel_in,T_fuel_out,burnup=0):
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


def fuel_restructuring(z,temp_fuel_out,temp_fuel_in,diam_fuel_out,diam_clad_out):
    """
    HP CONS:
    - using fv but neglecting beneficial increase of k_Fuel in col. region...
    - no equiaxial region designed

    GET:
    - new diameter of fuel by only ONE iteration (like hot geo iteration)
    - old diamtere of fuel
    - radius of columnar region
    - radius of void region
    - temperature @ void, hence new maximum: got using temp_fuel_max + void factor
    - old array of temperatures
    """
    radius_clmn = get_R_from_temp(z,diam_fuel_out,fuel_temp_clmn,temp_fuel_in,temp_fuel_out)
    radius_void = radius_void_get(radius_clmn,poro_asf)

    old = temp_fuel_max(z, diam_clad_out, diam_fuel_out, clad_thickness_0)
    if radius_clmn != 0:
        temp_void = temp_fuel_max(z,diam_clad_out,diam_fuel_out,clad_thickness_0,void_factor(radius_void,diam_fuel_out/2))
    else:
        temp_void = old

    ## hot geo - only one "iteration"
    old_diam = diam_fuel_out
    diam_fuel_out = diameter_th_exp_fuel(z, fuel_d_outer, temp_void, temp_fuel_out)

    return diam_fuel_out, old_diam, radius_clmn, radius_void, temp_void, old




#### **************************************************** GAS IN GAP ********************************************** ####

def gap_vol_cold():
    """
    Used to calculate initial moles of He in gap (as you'll see later, we'll also take into account the plenm extra vol)
    """
    vol_1 = pin_column_height * 0.25*pi*clad_d_outer**2
    vol_2 = pin_column_height * 0.25*pi*fuel_d_outer**2 # HP coincidenza altezze
    return vol_1 - vol_2


def gap_vol_hot(fuel_d_out_array,clad_d_out_array):
## check se corretto, inoltre tenere conto (In cold) dello spazio extra rihiesto per ex lungo asse
    num = len(fuel_d_out_array)
    unit = pin_column_height/num    # unit of length (discretization)
    vol = 0     # init volume
    clad_diameter_in = clad_d_out_array - 2*clad_thickness_0

    for i in range(0,num):
        # condition to avoid summing up a negative volume
        if clad_diameter_in[i] > fuel_d_out_array[i]:
            vol = vol + unit * 0.25 * pi * (clad_diameter_in[i] ** 2 - fuel_d_out_array[i] ** 2)
    return np.sum(vol)


def fg_prod(burnup):
    """
    Burn up in GWd/ton!!

    HP: FGR = 100 % (conservative/worst case and also easier)
    """
    N_av = 6.022e23 # atoms/mol
    y_xe = 0.27
    y_kr = 0.03
    y_he = 0.022

    E_fiss = 200e6 * 1.6e-19 # J
    q_third = power_lin_max/ (pi*0.25*fuel_d_outer**2) # W/m^3
    vol = pin_column_height * (pi*0.25*fuel_d_outer**2) # m^3
    F_prime = q_third / E_fiss # fiss/s/m^3

    mass = fuel_density * vol/1000 # ton
    q = power_lin_max * pin_column_height / 1e9 # GW
    time = burnup * 86400 * mass/ q # seconds

    mol_xe = ( y_xe * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    mol_kr = ( y_kr * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    mol_he = ( y_he * F_prime * vol ) / N_av  # atoms/s * mol/atoms = mol/s
    return mol_xe*time, mol_kr*time, mol_he*time


def pressure_gap_calc(volume_gap,temp_clad_array,burnup,plenum_vol=0,plenum_clad_d_in=0,temp_plenum=0):
    """
    Partial pressure computing and then summing them up
    HPs:
    - total release of gas (and also considering swelling, just see other section of the code)
    - temp plenum == coolant one
    - equation of perfect gas

    OSS:
    - temp_clad_array is matrix of temp along r and z respectively (rows and columns)
    - extra length in a
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
    if plenum_clad_d_in != 0: # condition to simply avoid zero at denominator :D
        plenum_length = plenum_vol / (0.25 * pi * plenum_clad_d_in ** 2) # m

    press_xe = (mol_xe * R * mean_temp) / vol
    press_kr = (mol_kr * R * mean_temp) / vol
    press_he = (mol_he * R * mean_temp) / vol
    total_pressure = press_he + press_xe + press_kr
    return total_pressure, plenum_length  # Pa

#temp = pressure_gap_calc(gap_vol_cold(), 800)
#print(f"{np.round(temp/1e6,2)} MPa")




#### ************************************************************************************************************** ####
#### ******************************************* PRINCIPAL FUNCTIONS ********************************************** ####
#### ************************************************************************************************************** ####

def thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0,bup):
    """
    Use it to get principal points temperatures: if used once, then you'll get cold geo temps, hence you should embed it
    in an iterative cycle properly set up in order to implement and get hot geometry temps.
    Also giving delta gap measures
    """
    output = np.zeros(5) # 0: temp coolant - 1: temp clad out - 2: temp clad in - 3: temp fuel out - 4: temp fuel in

    # cold geometry computing - to be used as first iteration data
    output[0] = temp_coolant(z)
    output[1] = temp_cladding_outer(z, clad_d_out_0)
    output[2] = temp_cladding_inner(z, clad_d_out_0, clad_thick_0)
    output[3], d_gap = temp_fuel_outer(z, clad_d_out_0, fuel_d_out_0, clad_thick_0, burnup=bup)
    output[4] = temp_fuel_max(z, clad_d_out_0, fuel_d_out_0, clad_thick_0, burnup=bup)
    return output, d_gap


def hot_geometry_general(z, clad_d_out_0, fuel_d_out_0, clad_thick_0,bup=0,print_status=True,RestructOn=True):
    """
    Iterative calculation in order to get temperatures, gap, final fuel out and clad out diameters, other properties
    (such as htc, net area, adim. numbers and so on...) in a hot geometry model


    NOTE:
    - Also accounting for fuel and cladding swelling
    - FOR RESTRUCTURING USEFUL UNTIL FUEL OUTER!
    - NEW! Implementation of burnup (our "time") eventually available, although by default set at zero (see above)
    - set RestructOn to take into account restructuring effects (however not computed if under a centain threshold)

    :param z: in m
    :param clad_d_out_0: cold geo diameter, in m
    :param fuel_d_out_0: cold geo diameter, in m
    :param clad_thick_0: supposed cladding thickness to be used, in m
    :param bup: Burnup in GWd/ton
    :return: cold temp, hot temp, minimun gap along the pin, generic properties @ hot geo, final clad outer diam, final fuel outer diam
    """
    tol = 10e-3 # tolerance set for following iteration phase (see later in this def)

    # 0: temp coolant - 1: temp clad out - 2: temp clad in - 3: temp fuel out - 4: temp fuel in
    temp_array,old_gap = thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0, bup)
    old_temp = temp_array.copy() # to give as output the cold geo temps too


    # variation of external diameters due to swelling
    # for cladding increase of thickness! (not more than about 4 %)
    fuel_d_outer = swelling_fuel(bup, fuel_d_out_0)
    clad_d_outer = swelling_clad(z, bup, temp_array[2], clad_d_out_0)
    old_clad_thick_0 = clad_thick_0
    clad_thick_0 += (clad_d_outer - clad_d_out_0)/2 # new thickness value computed

    ## CONSOLE PRINT SET UP here
    if print_status:
        print(f"** COMPLETED AT {np.round(100 * z / 0.85, 2)}% (POSITION: {np.round(z, 2)} m) **")
        print(f"GEO VARIATION DUE TO SWELLING:\nFuel diameter outer @ cold: {fuel_d_out_0} m --> Now: {fuel_d_outer} m")
        print(f"Cladding d outer @ cold: {clad_d_out_0} m --> Now: {round(float(clad_d_outer),8)} m")
        print(f"Clad thickness: {old_clad_thick_0} m --> Now: {round(float(clad_thick_0),8)} m")


    # START OF ITERATIVE COMPUTING OF TEMPERATURES #
    while True:
        prec_temp_array = temp_array.copy() # to be used to evaluate when exiting from the while below...

        # always consider expantion w.r.t. initial geometry (clad_d_outer, fuel_d_outer variables)
        clad_d_out_0 = diameter_th_exp_cladding(clad_d_outer, (prec_temp_array[2]))  # with temp clad inner HP CONS
        fuel_d_out_0 = diameter_th_exp_fuel(z, fuel_d_outer, prec_temp_array[4], prec_temp_array[3])

        # hot geo computing, starting by new geometries above
        temp_array, delta_gap = thermal_computing(z,clad_d_out_0, fuel_d_out_0, clad_thick_0,bup)

        if np.abs(prec_temp_array[4] - temp_array[4]) < tol:  # va bene così (?)
            break
    # END OF ITERATIVE COMPUTING OF TEMPERATURES #


    # preparing other stuff to be exported by the function
    yy_htc_loc, yy_adim_num_cool, yy_cool_loc_prop = heat_transfer_coefficient(temp_array[0], clad_d_out_0)  # new
    other = np.array(list([yy_htc_loc]) + list(yy_adim_num_cool) + list(yy_cool_loc_prop))

    # restructuring effects taken into account
    rst_threshold = 0.01 # GWd/ton
    if RestructOn and bup >= rst_threshold:
        print(f"RESTRUCTURING CONDITIONS SATISFIED!")
        fuel_d_out_0, _, _, temp_array[4], _,_ = fuel_restructuring(z,temp_array[3],temp_array[4],fuel_d_out_0,clad_d_out_0)

    ## CONSOLE PRINT SET UP here
    if print_status: # optional "progress bar" print (see input boolean)
        print(f"Fuel temp inner: HOT:{np.round(temp_array[4], 2)} K, COLD:{np.round(old_temp[4], 2)} K\nOld gap: {round(float(old_gap*1e6),4)} um --> New gap:{round(float(delta_gap*1e6),2)} um"
              f" ({round(float(100 * delta_gap / initial_delta_gap),2)}%)\n")

    return old_temp, temp_array, delta_gap, other, clad_d_out_0, fuel_d_out_0
