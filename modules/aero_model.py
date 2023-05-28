import numpy as np

def calc_C_thrust_x(phi, sigma, alpha, J):
    sin_alpha = np.sin(np.deg2rad(alpha))
    
    return phi * (phi / sigma * sin_alpha - J)

def calc_C_thrust_y(phi, sigma, alpha, J):
    cos_alpha = np.cos(np.deg2rad(alpha))
    
    return phi**2 / sigma * cos_alpha

def calc_C_power(phi, psi):
    return phi * psi

def calc_psi(phi, Cpr, J, eta):
    a = phi**2 * (1 - Cpr) - J**2

    return a / (2 * eta)


def calc_C_drag_x(alpha, L_D, hub_tip):
    # Geometry
    D2_A3 = 1 / (4 * np.pi) / (1 - hub_tip**2)

    # Drag
    f = 1.28

    sin_alpha = np.sin(np.deg2rad(alpha))
    cos_alpha = np.cos(np.deg2rad(alpha))

    return f * D2_A3 * (0.25 * np.pi * sin_alpha + L_D * cos_alpha)

def calc_C_drag_y(alpha):
    return 2 * np.pi * np.sin(np.deg2rad(alpha))