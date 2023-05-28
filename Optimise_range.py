from modules import geometry
from modules import aero_model
from modules import esdu_75026 as esdu
import scipy.optimize
from modules import propulsor_mass
import numpy as np
import matplotlib.pyplot as plt







# def calc_phi(J, sigma, alpha, L_D):

#     C_Tx_func = lambda phi: aero_model.calc_C_thrust_x(phi, sigma, alpha, J)
#     C_Dx = aero_model.calc_C_drag_x(alpha, L_D)

#     func = lambda phi: C_Tx_func(phi) - C_Dx

#     scipy.optimize.root(func,)




# def calc_range(J, C_power, M_bat):
#     return J / C_power * M_bat

# def calc_M_prop(sigma, geom):

#     # Diffuser params
#     Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)
#     angle = esdu.calc_phi(sigma, Lh)

#     # Update geometry
#     geom.phi_i = angle
#     geom.phi_o = angle
#     geom.Lh_diffuser = Lh

#     # Calc mass
    


# Choose initial parameters
M_payload = 12

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
    hub_solidity        =   0.2
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
# sigma = 1.2
eta = 0.9
# alpha = 10
# J = 0.5
# FrU = 5
phi = 0.7
hub_tip = geom.rr_hub_tip
L_D = geom.Lh_total / (2 * geom.Rh_cas)

def calc_s(sigma, alpha, FrU):

    # Update geometry for starting sigma
    Cpr, Lh = esdu.optimum_Cpr2_symmetrical(sigma)
    angle = esdu.calc_phi(sigma, Lh)

    Cpr = 1 - 1/sigma**2
    Lh *= 0.2

    geom.phi_i = angle
    geom.phi_o = angle
    geom.Lh_diffuser = Lh

    # Update mass model
    prop_mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)

    # Find J from x equation
    C_Tx_func = lambda J: aero_model.calc_C_thrust_x(phi, sigma, alpha, J)
    C_Dx = aero_model.calc_C_drag_x(alpha, L_D, hub_tip)
    f = lambda J: C_Tx_func(J) - C_Dx * J**2/2

    sol = scipy.optimize.root_scalar(f, bracket=(0, 20))
    # sol = scipy.optimize.root_scalar(f, x0=1, method='newton', fprime=fprime)
    J = sol.root

    # Find other variables
    psi = aero_model.calc_psi(phi, Cpr, J, eta)
    C_Wdotx = max(aero_model.calc_C_power(phi, psi), 0.06)
    C_Ty = aero_model.calc_C_thrust_y(phi, sigma, alpha, J)
    C_Dy = aero_model.calc_C_drag_y(alpha)

    M_total = (C_Ty - C_Dy * J**2 / 2) * FrU**2

    M_prop = prop_mass_model.calc_M_total()

    s = J / C_Wdotx * (M_total - M_prop - M_payload)

    return s


min_func = lambda x: -calc_s(*x)

# sigma, alpha, FrU
x0 = [1.2, 10, 5]
bounds = [(1, 3), (0, 50), (0, 60)]

# res = scipy.optimize.minimize(min_func, x0, method="Nelder-Mead", bounds=bounds)
res = scipy.optimize.minimize(min_func, x0, bounds=bounds)

print(-res.fun)
print(res.x)

sigma, alpha, FrU = res.x

# Find J from x equation
C_Tx_func = lambda J: aero_model.calc_C_thrust_x(phi, sigma, alpha, J)
C_Dx = aero_model.calc_C_drag_x(alpha, L_D, hub_tip)
f = lambda J: C_Tx_func(J) - C_Dx * J**2/2

sol = scipy.optimize.root_scalar(f, bracket=(0, 20))
J = sol.root

print(f"J = {J}")


# sigma = np.linspace(1, 2, 50)
# FrU = 60
# alpha = 10
# s = np.zeros_like(sigma)
# for i, sig in enumerate(sigma):
#     s[i] = calc_s(sig, alpha, FrU)

# plt.plot(sigma, s)
# plt.show()