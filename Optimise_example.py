from modules import geometry
from modules import aero_model
from modules import esdu_75026 as esdu
from modules import optimise
from modules import propulsor_mass
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plt.rcParams.update({"font.size": 20})



# ------------------------------------------------------------
#                Optimiser fixed parameters
# ------------------------------------------------------------

# Mass parameters
rho_air = 1.2
rho_nylon = 1150
rho_cu = 8960
rho_al = 2700
rho_resin = 1200

mass_params = propulsor_mass.NonDimensionalMassParams(
    cowl_density        =   rho_nylon / rho_air,
    hub_density         =   rho_cu / rho_air,
    rotor_density       =   rho_al / rho_air,
    stator_density      =   rho_resin / rho_air,
    cowl_solidity       =   0.3,
    hub_solidity        =   0.8
)

# Base geometry
geom = geometry.NonDimensionalGeometry(
    Lh_diffuser=0,    
    Lh_rotor=0.3,     
    Lh_intake=0.2,   
    rr_hub_tip=0.3, 
    th_cowl=0.3,
    tc_rotor=0.1,
    tc_stator=0.1,
    phi_i=0,
    phi_o=0,
    cL_rotor=0.9,
    cL_stator=0.6,
    delta_oh=0.02,
    N_rotor=6,
    N_stator=7
)

# Initial values/fixed values

# Diffuser Cpr and L/h func
def diffuser_Cpr_Lh_func(sigma):
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)

    # Cpr = 1 * (1 - 1/sigma**2)
    # Lh *= 0.3

    Lh = max(Lh, 0.2)

    return Cpr, Lh


# Dimensional parameters
M_payload_tot = 4
D = 0.12

FrU = 50
g = 9.8
Pa = 1e5
Ta = 15+273
R = 287
eta = 0.9
phi = 0.7
rhoE = 0.5e6

# Calculate dimensional parameters
M_payload = M_payload_tot / 4

Rc = D / 2
Rh = geom.rr_hub_tip * Rc
Rm = (Rc + Rh) / 2
A3 = np.pi * (Rc**2 - Rh**2)

rpm = np.sqrt(g * D) / (2 * np.pi * Rm) * FrU * 60
U = np.sqrt(g * D) * FrU

rho = Pa / (R * Ta)

M_payload_nd = M_payload / (rho * A3 * D)


# Optimiser
# sigma, alpha, FrU
x0 = [1.2, 20, 30]
bounds = [(1, 3), (0, 50), (0, FrU)]

optimise.diffuser_Cpr_Lh_func = diffuser_Cpr_Lh_func
res = optimise.optimise(x0, bounds, geom, mass_params, M_payload_nd, phi, eta)

names = ["sigma", "alpha", "FrU", "J"]


print()
for name in names:
    val = res.__dict__[name]
    print(f"{name} = {val}")
print()

# Post process
V = res.J * U
s = res.s * rhoE * D / U**2
M_bat = res.M_bat * (rho * A3 * D)
M_bat_tot = 4 * M_bat
M_prop = res.M_prop * (rho * A3 * D)
power = res.C_Wdotx * (rho * A3 * U**3)
E_bat_tot = M_bat_tot * rhoE / 3600

names = ["M_payload_tot", "rpm", "V", "s", "M_bat_tot", "M_prop", "Wdotx", "E_bat_tot"]
values = [M_payload_tot, rpm, V, s, M_bat_tot, M_prop, power, E_bat_tot]
units = ["kg", "rpm", "m/s", "m", "kg", "kg", "W", "Wh"]

print()
for name, val, unit in zip(names, values, units):
    print(f"{name} = {val} {unit}")
print()

fig, ax = plt.subplots(figsize=(10, 3))

hr = Rc - Rh
xr_hub = geom.calc_hub_line() * hr * 1000
xr_cas = geom.calc_cas_line() * hr * 1000
xr_cowl = geom.calc_cowl_line() * hr * 1000
curves = [xr_hub, xr_cas, xr_cowl]

for c in curves:
    ax.plot(c[1, :], -c[0, :], '-r')
    ax.plot(-c[1, :], -c[0, :], '-r')

ax.set_aspect('equal')
ax.grid(True, which='both')
ax.set_xlabel("$ r $ (mm)")
ax.set_ylabel("$ x $ (mm)", rotation=0, labelpad=20)
ax.yaxis.set_minor_locator(mticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(mticker.MultipleLocator(5))
ax.xaxis.set_major_formatter(lambda x, pos: str(int(abs(x))))
fig.set_constrained_layout(True)

fig.savefig("Figures/Optimise_example_geom.pdf")
plt.show()

