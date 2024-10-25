from general_functions import *
from material_properties import *
from design_specifications import *


#### TEMP PROFILE OF COOLANT ALONG Z AXIS ####
def temp_coolant(z):
    """
    :param z: position choice
    :return: Temp of coolant at point T
    """
    # energy balance: f*int(Cp(T), Tin, T) - int(q', bottom_pin, z) / mass_flow = 0 ... Tout our goal.
    # occorre metodo che minimizzi la funzione res affinche si trovi una temperatura che rispetti il bilancio energetico, in quanto Cp stesso dipende da temperatura!! x primo membro
    T_0 = 800 + 273.15 # K
    eqz_2 = f_pitch*integral_power_lin_distr(z) / cool_mass_flow  # valore num seconda parte eqz !! NB F PITCH=1/2 !!
    eqz_1 = sy.integrate(cool_spec_heat,temp)
    res = eqz_1 - eqz_1.subs(temp,cool_temp_inlet) - eqz_2 # primo è integtale, secondo intehrale con sostituzione (quindi risolto) e terzo eqz_2

    T_coolant_at_z = equation_temp_solver(res, T_0)
    return T_coolant_at_z


#### CORRELATION FUNCTIONS ####
# in un punto preciso (esso può variare lungo z ????????) in base a temp in K; poi aggiungere dip da hydr diameter per HOT GEOMETRY
# {NB: dubbi su htc, occorre verifica h e Pr}
def heat_transfer_coeff_local(temperature,clad_diam_out):
    """
    :param clad_diam_out: in m - Varying parameter! (hot geometry)
    :param temperature: in [K]
    :return:
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
    return h, output_numbers


#### TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS ####
# un giorno modifica pure cladding outer diam (th exp)
def temp_cladding_outer(z,clad_d_out):
    """
    :param z: in m
    :param clad_d_out: in m - Varying parameter! (hot geometry)
    :return: in K
    """
    sec_power_at_z = power_lin_distribution(z) / (pi*clad_d_out)
    temp_coolant_at_z = temp_coolant(z)
    htc,other = heat_transfer_coeff_local(temp_coolant_at_z,clad_d_out) # other contiene numeri adim, che non servono
    temp_clad_out = temp_coolant_at_z + sec_power_at_z/htc
    return temp_clad_out


### here unknown thickness!!! inoltre hot geometry -> variazione D_clad_outer; HP cons delta thickness cost...
def temp_cladding_inner(z,clad_d_out,clad_thick):
    """
    :param clad_d_out: Varying parameter! (hot geometry)
    :param z: in [m]
    :param clad_thick: in [m] - UNKNOWN
    :return: in K
    """
    temp_ci_guess = (600+273.15)

    eqz_1 = temp - temp_cladding_outer(z,clad_d_out) # variable: temp presa da mat_properties...
    eqz_2 = power_lin_distribution(z) * clad_thick / ( pi * (clad_d_outer-2*clad_thick) * clad_thermal_cond ) # clad th cond dep. on temp too!
    res = eqz_1 - eqz_2

    output = equation_temp_solver(res, temp_ci_guess)
    return output


# NB variazione diam ext del fuel (hot geometry)
def temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th):
    """
    DELTA GAP UNKNOWN, it depends on Fuel out and Clad in, which depends itself on Thickness (CONST)
        and Clad out (fixed by project but expanding...)
    :param clad_d_out: Varying parameter! (hot geometry)
    :param z: in m
    :param fuel_diam_outer: in m - Varying parameter! (hot geometry)
    :param clad_th: in m - UNKNOWN
    :return: in K
    """
    temp_fuel_outer_guess = 1000 + 273.15
    delta_gap = (clad_d_out - 2 * clad_th - fuel_diam_outer) / 2
    delta_gap_eff = delta_gap + 10e-6 # m - 10e-6 He poiché effettivo...
    temp_clad_in = temp_cladding_inner(z,clad_d_out,clad_th)

    eqz_1 = temp - temp_clad_in
    eqz_2 = power_lin_distribution(z) * delta_gap_eff / ( pi * fuel_diam_outer * helium_thermal_cond )
    res = eqz_1 - eqz_2

    out = equation_temp_solver(res, temp_fuel_outer_guess)
    return out


### da espandere per bene (void factor, zone restructuring, pu redistri... è una bozza al momento!!)
def temp_fuel_inner(z,clad_d_out,fuel_diam_outer,clad_th):
    """
    :param z: in m
    :param clad_d_out: - Varying parameter! (hot geometry)
    :param fuel_diam_outer: - Varying parameter! (hot geometry)
    :param clad_th: UNKNOWN
    :param delta_gap: - Varying parameter! (hot geometry) - UNKNOWN
    :return: in K
    """
    temp_fuel_inner_guess = 1200 + 273.15
    temp_fuel_out = temp_fuel_outer(z,clad_d_out,fuel_diam_outer,clad_th)

    k_fuel = fuel_thermal_cond.subs(x_om,2) # per ora questi valori
    k_fuel = k_fuel.subs(pu_conc,0.2)
    k_fuel = k_fuel.subs(por,0.12)

    eqz_1 = temp - temp_fuel_out
    eqz_2 = power_lin_distribution(z) / (4*pi*k_fuel)
    #print(k_fuel)

    res = eqz_1 - eqz_2
    out = equation_temp_solver(res, temp_fuel_inner_guess)
    return out


def temp_fuel_inner_radial(z,clad_d_out_fuel_diam_outer,clad_th):



    return None


#### HOT GEOMETRY FUNCTION ####
def diameter_th_exp_cladding(diam,temperature):
    #return diam + diam * clad_eps_th.subs(temp,temperature) * (temperature - temp_in)
    return diam + diam * clad_eps_th * (temperature - temp_in) ### RIVEDI QUESTA ROBA CHE FORSE MANCA IL DATO!

def diameter_th_exp_fuel(diam,temperature):
    return diam + diam * alfa_fuel * (temperature - temp_in)