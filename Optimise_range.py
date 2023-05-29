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

# Payload mass
M_payload = 10

# Diffuser Cpr and L/h func
def diffuser_Cpr_Lh_func(sigma):
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)

    Cpr = 1 * (1 - 1/sigma**2)
    Lh *= 0.3

    return Cpr, Lh


# sigma, alpha, FrU
x0 = [1.2, 20, 30]
bounds = [(1, 3), (0, 50), (0, 60)]

optimise.diffuser_Cpr_Lh_func = diffuser_Cpr_Lh_func
res = optimise.optimise(x0, bounds, geom, mass_params, M_payload, phi, eta)

for param, val in res.__dict__.items():
    if not isinstance(val, float):
        continue

    print(f"{param} = {val}")

print(f"M_payload/M_total (%) = {100*M_payload/res.M_total}")
print(f"M_bat/M_total (%) = {100*res.M_bat/res.M_total}")


fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(res.s_conv / res.s_conv[-1])
ax.plot(res.sigma_conv / res.sigma_conv[-1])
ax.plot(res.alpha_conv / res.alpha_conv[-1])
ax.plot(res.FrU_conv / res.FrU_conv[-1])
ax.grid()
ax.legend([r"$ \widetilde s $", r"$ \sigma $", r"$ \alpha $", r"$ Fr_U $"])
ax.set_xlabel("Number of iterations")
ax.set_ylabel(r"$ \dfrac{( \ )}{( \ )_{opt}} $", rotation=0, labelpad=20)
ax.set_title(r"$ M_{payload} = " + str(M_payload) + r" $")
fig.set_constrained_layout(True)

response = input("Save file? (y/n): ")
if response == "y":
    fig.savefig(f"Figures/Optimise_conv_Mpl{M_payload}_both.pdf")

plt.show()
