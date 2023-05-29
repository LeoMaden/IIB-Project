from modules import geometry
from modules import aero_model
from modules import esdu_75026 as esdu
from modules import optimise
from modules import propulsor_mass
import numpy as np
import matplotlib.pyplot as plt

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
    cowl_solidity       =   0.1,
    hub_solidity        =   0.25
)

# Base geometry
geom = geometry.NonDimensionalGeometry(
    Lh_diffuser=0.3,    
    Lh_rotor=0.2,     
    Lh_intake=0.1,   
    rr_hub_tip=0.15, 
    th_cowl=0.3,
    tc_rotor=0.1,
    tc_stator=0.1,
    phi_i=20,
    phi_o=20,
    cL_rotor=0.9,
    cL_stator=0.6,
    delta_oh=0.02,
    N_rotor=6,
    N_stator=7
)


# Initial values/fixed values
eta = 0.9
phi = 0.7

# Diffuser Cpr and L/h func
def diffuser_Cpr_Lh_func(sigma):
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)

    # Cpr = 1 * (1 - 1/sigma**2)
    # Lh *= 0.3

    return Cpr, Lh


# sigma, alpha, FrU
x0 = [1.2, 20, 30]
bounds = [(1, 3), (0, 50), (0, 60)]

optimise.diffuser_Cpr_Lh_func = diffuser_Cpr_Lh_func

# Payload mass
M_payload_arr = np.linspace(0, 1000, 100)
results = []

for i, M_payload in enumerate(M_payload_arr):
    res = optimise.optimise(x0, bounds, geom, mass_params, M_payload, phi, eta)
    results.append(res)


# Parameters
fig1, ax1 = plt.subplots(figsize=(10, 8))

params = ["s", "sigma", "alpha", "J", "psi"]
p_names = ["s", r"\sigma", r"\alpha", "J", r"\psi"]

for p in params:

    val0 = results[0].__dict__[p]

    val_arr = [res.__dict__[p]/val0 for res in results]

    ax1.plot(M_payload_arr, val_arr)

ax1.legend([f"$ {p} $" for p in p_names])
ax1.set_xlabel(r"$ \widetilde M_{payload} $")
ax1.grid('major')

fig1.set_constrained_layout(True)

# Mass
fig2, ax2 = plt.subplots(figsize=(10, 8))

M_total_arr = np.array([res.M_total for res in results])
M_prop_arr = np.array([res.M_prop for res in results])
M_bat_arr = np.array([res.M_bat for res in results])

M_payload_rel = M_payload_arr #/ M_total_arr
M_prop_rel = M_prop_arr #/ M_total_arr
M_bat_rel = M_bat_arr #/ M_total_arr

stack = np.vstack([M_bat_rel, M_prop_rel, M_payload_rel])

ax2.stackplot(M_payload_arr, stack)

ax2.legend([r"$ \widetilde M_{bat} $", r"$ \widetilde M_{prop} $", r"$ \widetilde M_{payload} $"], loc='lower right')
ax2.set_xlabel(r"$ \widetilde M_{payload} $")
ax2.grid('major')

fig2.set_constrained_layout(True)


fig1.savefig("Figures/Optimise_vs_payload_params_ESDU.pdf")
fig2.savefig("Figures/Optimise_vs_payload_mass_ESDU.pdf")
# fig1.savefig("Figures/Optimise_vs_payload_params_both.pdf")
# fig2.savefig("Figures/Optimise_vs_payload_mass_both.pdf")

plt.show()


