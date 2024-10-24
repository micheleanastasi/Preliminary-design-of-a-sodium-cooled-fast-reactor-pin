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

    T_coolant_at_z = iterative_solver(res,T_0)
    return T_coolant_at_z


def plot_temp_coolant():

    xx = np.linspace(pin_bottom_pos, pin_top_pos, 100)
    yy = np.zeros(100)

    for i in range(0,len(xx)):
        yy[i] = temp_coolant(xx[i])
    plt.figure("Coolant temperature profile")
    plt.plot(xx,yy,label = 'Temperature profile')
    plt.xlabel("Fuel pin position in [m]")
    plt.ylabel("Temperature in [K]")
    plt.title("Coolant temperature profile")
    plt.legend()
    plt.grid()
    plt.show()
    return None


#### CORRELATION FUNCTIONS

# in un punto preciso (esso può variare lungo z ????????) in base a temp in K; poi aggiungere dip da hydr diameter per HOT GEOMETRY
def heat_transfer_coeff_local(temperature):

    temp_fahr = (temperature - 273.15) * 9 / 5 + 32 # conversione kelvin to fahrenheit
    density = cool_density__fahr.subs(temp_f,temp_fahr) # kg/m^3
    dyn_viscosity = cool_dyn_viscosity.subs(temp,temperature) # Pa*s
    spec_heat = cool_spec_heat.subs(temp,temperature) # J/Kg/K
    th_cond = cool_th_cond.subs(temp,temperature) # J/m/K

    net_area = pin_pitch**2 * sqrt(3)/4 - 0.5*( pi*clad_d_outer**2/4 ) # m^2
    hydr_diameter = clad_d_outer * ( (2*sqrt(3)/pi)*(pin_pitch/clad_d_outer)**2 - 1 ) # m
    avg_velocity = cool_mass_flow/( net_area*density) # m/s

    Re = density*avg_velocity*hydr_diameter/dyn_viscosity
    Pr = spec_heat*dyn_viscosity/th_cond
    Pe = Pr*Re

    Nu = Nu_cool.subs(Pe_cool,Pe)

    h = Nu*th_cond/hydr_diameter
    output_numbers = array([Re,Pr,Pe,Nu])
    return h, output_numbers


#### TEMP PROFILE OF CLADDING (OUTER) ALONG Z AXIS ####
# un giorno modifica pure cladding outer diam (th exp)
def temp_cladding_outer(temp_coolant_at_z,lin_power_at_z,correlation):

    sec_power_at_z = lin_power_at_z / (pi*clad_d_outer)
    temp_clad_out = temp_coolant_at_z + sec_power_at_z/correlation

    return temp_clad_out

def plot_temp_cladding_outer():
    xx = np.linspace(pin_bottom_pos, pin_top_pos, 100)
    yy = np.zeros(100)
    for i in range(0,len(xx)):
        temp_at_z = temp_coolant(xx[i])
        correlation = heat_transfer_coeff_local( temp_at_z )[0]
        power_lin = power_lin_distribution(xx[i])
        yy[i] = temp_cladding_outer(temp_at_z, power_lin, correlation)
    plt.figure("Coolant temperature profile")
    plt.plot(xx,yy,label = 'Temperature profile')
    #plt.plot(xx,yy_const,'.', label = 'Boiling temperature')
    plt.xlabel("Fuel pin position in [m]")
    plt.ylabel("Temperature in [K]")
    plt.title("Coolant temperature profile")
    plt.legend()
    plt.grid()
    plt.show()
    return None

### here unknown thickness!!! inoltre hot geometry, variazione D clad outer HP cons delta thickness cost...
def temp_cladding_inner(z,clad_thick):
    """
    :param z: in [m]
    :param clad_thick: in [m]
    :return:
    """
    temp_ci_guess = (600+273.15)
    temp_cool = temp_coolant(z)
    htc,other = heat_transfer_coeff_local(temp_cool)

    eqz_1 = temp - temp_cladding_outer( temp_cool, power_lin_distribution(z), htc ) # variable: temp
    eqz_2 = power_lin_distribution(z) * clad_thick / ( pi * (clad_d_outer-2*clad_thick) * clad_thermal_cond ) # clad th cond dep. on temp too!
    res = eqz_1 - eqz_2

    output = iterative_solver(res,temp_ci_guess)
    return output

"""
plotting(temp_coolant)
test = temp_coolant(0.85) - 273.15
print(f"T all'uscita del coolant in °C: {np.round(test,2)}")


h, numbers = heat_transfer_coeff_local(400 + 273.15)
print(h)
print(numbers)

plot_temp_cladding_outer()
temp_at_z = temp_coolant(0.85)
correlation = heat_transfer_coeff_local(temp_at_z)[0]
power_lin = power_lin_distribution(0.85)
res = temp_cladding_outer(temp_at_z, power_lin, correlation)
res -= 273.15

print(res)
"""
print(temp_cladding_inner(0.85,80e-6))