from modules import geometry
from modules import aero_model
from modules import esdu_75026 as esdu
import scipy.optimize
from modules import propulsor_mass
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class OptimiseResult:
    # Result
    s: float

    # Primary
    sigma: float
    alpha: float
    FrU: float

    # Secondary
    J: float
    C_Tx: float
    C_Ty: float
    C_Dx: float
    C_Dy: float
    psi: float
    C_Wdotx: float

    # Mass
    M_prop: float
    M_bat: float
    M_total: float

    # Diffuser
    Cpr: float
    Lh: float

    # Output from optimiser
    res: scipy.optimize.OptimizeResult

    # Convergence history
    s_conv: np.array
    sigma_conv: np.array
    alpha_conv: np.array
    FrU_conv: np.array
    



diffuser_Cpr_Lh_func = esdu.optimum_Cpr2_symmetrical



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

    return J, C_Tx, C_Ty, C_Dx, C_Dy, psi, C_Wdotx

def calc_new_geom(sigma, geom):

    # Update geometry for starting sigma
    Cpr, Lh = diffuser_Cpr_Lh_func(sigma)

    # Calc new geometry
    angle = esdu.calc_phi(sigma, Lh)
    geom.phi_i = angle
    geom.phi_o = angle
    geom.Lh_diffuser = Lh

    return Cpr, Lh

def calc_masses(C_Ty, C_Dy, J, FrU, M_payload, prop_mass_model):
    M_total = (C_Ty - C_Dy * J**2 / 2) * FrU**2
    M_prop = prop_mass_model.calc_M_total()
    M_bat = (M_total - M_prop - M_payload)

    return M_prop, M_bat, M_total

def calc_s(sigma, alpha, FrU, phi, eta, M_payload, geom, mass_params, L_D, hub_tip):

    # Update geometry
    Cpr, Lh = calc_new_geom(sigma, geom)

    # Update mass model
    prop_mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)

    # Calculate model
    J, C_Tx, C_Ty, C_Dx, C_Dy, psi, C_Wdotx = calc_aero(phi, sigma, alpha, Cpr, eta, L_D, hub_tip)

    # Calculate mass
    M_prop, M_bat, M_total = calc_masses(C_Ty, C_Dy, J, FrU, M_payload, prop_mass_model)

    # Calculate range
    s = J / C_Wdotx * M_bat

    return s



def _calc_s_conv(sigma, alpha, FrU, phi, eta, M_payload, geom, mass_params, L_D, hub_tip, s_conv, sigma_conv, alpha_conv, FrU_conv):
    s = calc_s(sigma, alpha, FrU, phi, eta, M_payload, geom, mass_params, L_D, hub_tip)

    s_conv.append(s)
    sigma_conv.append(sigma)
    alpha_conv.append(alpha)
    FrU_conv.append(FrU)

    return s




def optimise(x0, bounds, geom, mass_params, M_payload, phi, eta):

    s_conv = []
    sigma_conv = []
    alpha_conv = []
    FrU_conv = []
    
    hub_tip = geom.rr_hub_tip
    L_D = geom.Lh_total / (2 * geom.Rh_cas)

    s_func = lambda x: _calc_s_conv(*x, phi, eta, M_payload, geom, mass_params, L_D, hub_tip, s_conv, sigma_conv, alpha_conv, FrU_conv)
    obj_func = lambda x: -s_func(x)

    res = scipy.optimize.minimize(obj_func, x0, method='L-BFGS-B', bounds=bounds)

    sigma, alpha, FrU = res.x
    s_opt = -res.fun

    # Calculate secondary
    Cpr, Lh = calc_new_geom(sigma, geom)
    prop_mass_model = propulsor_mass.NonDimensionalMassModel(mass_params, geom)
    J, C_Tx, C_Ty, C_Dx, C_Dy, psi, C_Wdotx = calc_aero(phi, sigma, alpha, Cpr, eta, L_D, hub_tip)
    M_prop, M_bat, M_total = calc_masses(C_Ty, C_Dy, J, FrU, M_payload, prop_mass_model)

    out = OptimiseResult(
        s=s_opt,
        sigma=sigma,
        alpha=alpha,
        FrU=FrU,
        J=J,
        C_Tx=C_Tx,
        C_Ty=C_Ty,
        C_Dx=C_Dx,
        C_Dy=C_Dy,
        psi=psi,
        C_Wdotx=C_Wdotx,
        M_prop=M_prop,
        M_bat=M_bat,
        M_total=M_total,
        Cpr=Cpr,
        Lh=Lh,
        res=res,
        s_conv=np.array(s_conv),
        sigma_conv=np.array(sigma_conv),
        alpha_conv=np.array(alpha_conv),
        FrU_conv=np.array(FrU_conv)
    )

    return out