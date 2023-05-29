from modules import geometry
from modules import aero_model
from modules import esdu_75026 as esdu
import scipy.optimize
from modules import propulsor_mass
import numpy as np
import matplotlib.pyplot as plt



# ------------------------------------------------------------
#                Optimiser fixed parameters
# ------------------------------------------------------------



# Choose initial parameters
M_payload = 10

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



# ------------------------------------------------------------



def calc_aero(phi, sigma, alpha, Cpr, eta, L_D, hub_tip):
    vals = {}

    # Find J from x equation
    C_Tx_func = lambda J: aero_model.calc_C_thrust_x(phi, sigma, alpha, J)
    C_Dx = aero_model.calc_C_drag_x(alpha, L_D, hub_tip)
    f = lambda J: C_Tx_func(J) - C_Dx * J**2/2

    sol = scipy.optimize.root_scalar(f, bracket=(0, 20))
    J = sol.root

    # Find other variables
    C_Tx = C_Tx_func(J)
    psi = aero_model.calc_psi(phi, Cpr, J, eta)
    C_Wdotx = aero_model.calc_C_power(phi, psi)
    C_Ty = aero_model.calc_C_thrust_y(phi, sigma, alpha, J)
    C_Dy = aero_model.calc_C_drag_y(alpha)

    vals["J"] = J
    vals["C_Tx"] = C_Tx
    vals["C_Ty"] = C_Ty
    vals["C_Dx"] = C_Dx
    vals["C_Dy"] = C_Dy
    vals["psi"] = psi
    vals["C_Wdotx"] = C_Wdotx

    return vals

def calc_new_geom(sigma, geom):

    # Update geometry for starting sigma
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)

    # Manipulate Cpr or L/h
    Cpr = 1 * (1 - 1/sigma**2)
    Lh *= 0.3

    # Calc new geometry
    angle = esdu.calc_phi(sigma, Lh)
    geom.phi_i = angle
    geom.phi_o = angle
    geom.Lh_diffuser = Lh

    return Cpr, Lh

def calc_masses(C_Ty, C_Dy, J, FrU, prop_mass_model):
    M_total = (C_Ty - C_Dy * J**2 / 2) * FrU**2
    M_prop = prop_mass_model.calc_M_total()
    M_bat = (M_total - M_prop - M_payload)

    return M_prop, M_bat, M_total


hub_tip = geom.rr_hub_tip
L_D = geom.Lh_total / (2 * geom.Rh_cas)

sigma_conv = []
alpha_conv = []
FrU_conv = []
s_conv = []

def calc_s(sigma, alpha, FrU):

    sigma_conv.append(sigma)
    alpha_conv.append(alpha)
    FrU_conv.append(FrU)

    # Update geometry
    Cpr, Lh = calc_new_geom(sigma, geom)

    # Update mass model
    prop_mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)

    # Calculate model
    vals = calc_aero(phi, sigma, alpha, Cpr, eta, L_D, hub_tip)

    C_Ty = vals["C_Ty"]
    C_Dy = vals["C_Dy"]
    C_Wdotx = vals["C_Wdotx"]
    J = vals["J"]

    # Calculate mass
    M_prop, M_bat, M_total = calc_masses(C_Ty, C_Dy, J, FrU, prop_mass_model)

    # Calculate range
    s = J / C_Wdotx * M_bat

    s_conv.append(s)
    return s


min_func = lambda x: -calc_s(*x)

# sigma, alpha, FrU
x0 = [1.2, 10, 30]
bounds = [(1, 3), (0, 50), (0, 60)]

# res = scipy.optimize.minimize(min_func, x0, method="Nelder-Mead", bounds=bounds)
res = scipy.optimize.minimize(min_func, x0, method='L-BFGS-B', bounds=bounds)

sigma, alpha, FrU = res.x

print(res.message)
print("s = {:.1f}".format(-res.fun))
print(f"sigma = {sigma:.3g}, alpha = {alpha:.3g}, FrU = {FrU:.3g}")


# Print aero values at end
Cpr, Lh = calc_new_geom(sigma, geom)
prop_mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)
vals = calc_aero(phi, sigma, alpha, Cpr, eta, L_D, hub_tip)
M_prop, M_bat, M_total = calc_masses(vals["C_Ty"], vals["C_Dy"], vals["J"], FrU, prop_mass_model)

for key, val in vals.items():
    print(f"{key} = {val:.3g}")

print(f"Cpr = {Cpr:.3g}")
print(f"Ld/hr = {Lh:.3g}")
print(f"M_payload = {M_payload:.3g}")
print(f"M_prop = {M_prop:.3g}")
print(f"M_bat = {M_bat:.3g}")
print(f"M_total = {M_total:.3g}")

plt.plot(sigma_conv/sigma_conv[-1])
plt.plot(alpha_conv/alpha_conv[-1])
plt.plot(FrU_conv/FrU_conv[-1])
plt.plot(s_conv/s_conv[-1])
plt.grid()
plt.show()
