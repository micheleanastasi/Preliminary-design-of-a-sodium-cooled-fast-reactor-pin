"""
PRELIMINARY SIZING FOR REPORT
"""
import sympy as sym
import matplotlib.pyplot as plt

from functions.general_properties import *
from functions.mechanical_functions import *





##   Numerical data, SI units

pressure_gas_cold=0.1e6    #[Pa]
T_He_cold=20+273.15        #[K]
T_He_hot=1000              #[K], guess
R=8.3145                   #[J/mol*K]  - gas constant
lin_power=38.7e3           #[W/m]
core_length=0.85    #[m]
Pu_mass_ratio=0.29  #wt
time_irradiation=360*24*60*60   #[s]
energy_fission=200*1e6*1.6e-19  #[J]
N_A=6.022e23
neutron_flux = 6.1e15  # n cm^-2 s^-1 (peak neutron flux)

fission_gas_release=1 # assuming worst case scenario- all gas produced are released

y_xe = 0.27    #yield
y_kr = 0.03
y_he = 0.022

clad_thickness=0.53e-3 #[mm]

#   design limits

max_pressure_plenum=5e6 #[Pa]
max_clad_temperature=650+273.15 #[K]




##  Burnup

fuel_volume=np.pi*fuel_d_outer**2/4*core_length
fuel_mass=fuel_volume*fuel_density
burnup_MOX=lin_power*1e-9*core_length*(time_irradiation/60/60/24)/(fuel_mass*1e-3)
print(f'Burnup [GWd/ton_MOX] : {round(burnup_MOX,3)}')

molar_mass_MOX=2*16+Pu_mass_ratio*244.06+(1-Pu_mass_ratio)*238.0289 #[g/mol]

fuel_mass_HM=fuel_mass/molar_mass_MOX*(Pu_mass_ratio*244.06+(1-Pu_mass_ratio)*238.0289)   #   heavy metals only
burnup_HM=lin_power*1e-9*core_length*(time_irradiation/60/60/24)/(fuel_mass_HM*1e-3)
print(f'Burnup [GWd/ton_HM] : {round(burnup_HM,3)}')




#   Embrittlement

def yield_strength_cladding(T):
    if T < 600:
        return 555.5e6 - 0.25 * T
    elif T <= 1000:
        return 405.5e6 - 0.775 * (T - 600)
    else:
        return 345.5e6 - 0.25 * T

stress_yield=yield_strength_cladding(max_clad_temperature-273.15)



# 1. Helium embrittlement of cladding

helium_content_moles = helium_content(burnup_HM)   #He content in cladding

stress_yield_embrittled = yield_stress_embrittlement(stress_yield, helium_content_moles)

molar_mass_MOX=2*16+Pu_mass_ratio*244.06+(1-Pu_mass_ratio)*238.0289 #[g/mol]

#   assuming 100% released
fission_gas_moles=power_lin_avg*core_length*time_irradiation/energy_fission/N_A*(y_he+y_xe+y_kr)*fission_gas_release
pressure_He_hot=pressure_gas_cold*T_He_hot/T_He_cold

Sym_volume=sym.Symbol('V')
sym_moles=sym.Symbol('n')

volume_plenum=sym.Symbol('V')
He_moles=volume_plenum*pressure_gas_cold/R/T_He_cold




#   pressurization of the fuel rod: unknown volume is filled with He until p=0.1 [MPa], then sealed
#   still assuming cold geometry, as fission gas are released pressure builds up while V is constant
#   neglecting other effects: swelling, gap reduction

pressure_plenum=(He_moles+fission_gas_moles)*R*T_He_hot/volume_plenum
num_pressure_plenum=sym.lambdify(volume_plenum,pressure_plenum)

plt.plot(np.linspace(0,2,1000),num_pressure_plenum(np.linspace(0,2,1000))*1e-6)
plt.xlabel('plenum volume [m^3]')
plt.ylabel('pressure [MPa]')
plt.title('Plenum pressure vs plenum volume')
plt.grid()
plt.show()

min_volume_plenum=sym.solve(sym.Eq(pressure_plenum,max_pressure_plenum))[0]
print(f'Minumum volume (gap+plenum) : {min_volume_plenum} [m^3]')



#   cladding thickness - Mariotte solution - Tresca criterion

max_pressure=max_pressure_plenum+25e6   # assuming contact pressure of about 25 MPa
min_thickness=max_pressure*clad_d_outer/2/(2/3*stress_yield_embrittled+0.5*max_pressure_plenum)
print(f'minimum cladding thickness : {min_thickness*1e3} [mm]')

print(f"Ratio clad d outer / thickness: {clad_d_outer/clad_thickness}")


clad_d_inner=clad_d_outer-clad_thickness*2

#   worst case is the gap shrinkens to near 0 - compute plenum volume accordingly

min_plenum_height=min_volume_plenum/(np.pi*(clad_d_inner)**2/4)   #[m]

print(f'Minimum plenum height (before swelling): {min_plenum_height} [m]')

#   void swelling

fluence=neutron_flux*time_irradiation

def void_swelling(T, phi):
    term1 = 1.5e-4
    term2 = np.exp(-2500 / (T - 450))
    term3 = phi**0.16
    return term1 * term2 * term3

void_swelling_cladding = void_swelling(max_clad_temperature, fluence)
void_swelling_cladding_lin=void_swelling_cladding*clad_d_outer/2/3
print(f"Void swelling of cladding: {void_swelling_cladding:.6f} (fractional increase in volume)")

##  Preliminary sizing

# thickness: since the yielding stress is very low, we will size it according to the needs of the temperature profile, then verify a posteriori
# plenum height : use operating T, BEFORE void swelling, as swelling increases the inner radius

volume_safety_margin=0   #   adjustable

min_volume_plenum_wmargin=min_volume_plenum*(1+volume_safety_margin)
min_plenum_height_wmargin=min_volume_plenum_wmargin/(np.pi*(clad_d_inner)**2/4)
print(f'Minimum plenum height : {round(min_plenum_height_wmargin,4)} [m]')

gap_cold=clad_d_inner/2-fuel_d_outer/2
print(f'Gap size (cold) : {gap_cold*1e3} [mm]')