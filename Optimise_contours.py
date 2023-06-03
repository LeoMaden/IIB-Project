from modules import optimise
from modules import geometry
from modules import esdu_75026 as esdu
from modules import propulsor_mass
import matplotlib.pyplot as plt
import numpy as np

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
    
hub_tip = geom.rr_hub_tip
L_D = geom.Lh_total / (2 * geom.Rh_cas)


# Initial values/fixed values
eta = 0.9
phi = 0.7
FrU = 60

# Payload mass
M_payload = 400

# Diffuser Cpr and L/h func
def diffuser_Cpr_Lh_func(sigma):
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)

    Cpr = 1 * (1 - 1/sigma**2)
    Lh *= 0.3

    return Cpr, Lh

optimise.diffuser_Cpr_Lh_func = diffuser_Cpr_Lh_func


# ------------------------------------------------------------


n = 50

alpha = np.linspace(0, 40, n)
sigma = np.linspace(1, 2, n)
Alpha, Sigma = np.meshgrid(alpha, sigma)

s = np.zeros((n, n))


for i in range(n):
    for j in range(n):
        alpha_ij = Alpha[i, j]
        sigma_ij = Sigma[i, j]
        s[i, j] = optimise.calc_s(sigma_ij, alpha_ij, FrU, phi, eta, M_payload, geom, mass_params, L_D, hub_tip)

fig, ax = plt.subplots(figsize=(10, 8))

max_level = np.ceil(np.max(s) / 50) * 50
levels = np.arange(0, max_level, 50)

cf = ax.contourf(Alpha, Sigma, s, cmap='plasma', levels=levels)
cbar = fig.colorbar(cf)
cbar.set_label("$ \widetilde s $", rotation=0, labelpad=20)
# cbar.

ax.set_xlabel(r"$ \alpha $ (deg)")
ax.set_ylabel(r"$ \sigma $", rotation=0, labelpad=20)
# ax.set_title(r"$ \widetilde M_{payload} = " + str(M_payload) + r" $")

fig.savefig(f"Figures/Optimise_contours_Mpl{M_payload}_Ypi1.pdf")

plt.show()
quit()


# ------------------------------------------------------------


# sigma, alpha, FrU
x0 = [1.2, 20, 30]
bounds = [(1, 3), (0, 50), (0, 60)]

# Payload mass
M_payload_arr = np.linspace(0, 300, 301)
results = []

for i, M_payload in enumerate(M_payload_arr):
    res = optimise.optimise(x0, bounds, geom, mass_params, M_payload, phi, eta)
    results.append(res)

# Parameters
fig2, ax2 = plt.subplots(figsize=(10, 8))

params = ["s", "sigma", "alpha"]
p_names = ["s", r"\sigma", r"\alpha"]

for p in params:

    val0 = results[0].__dict__[p]

    val_arr = [res.__dict__[p]/val0 for res in results]

    ax2.plot(M_payload_arr, val_arr)

ax2.legend([f"$ {p} $" for p in p_names])
ax2.set_xlabel(r"$ \widetilde M_{payload} $")
ax2.set_ylabel(r"$ ( \ ) \div ( \ ) $ at zero payload")
ax2.grid('major')

fig2.set_constrained_layout(True)

fig2.savefig("Figures/Optimise_contours_payload_variation.pdf")

plt.show()