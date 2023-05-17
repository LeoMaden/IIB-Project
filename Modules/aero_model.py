import numpy as np

# class AeroModel:
#     pass


def calc_C_thrust_x(phi, sigma, alpha, J, P1):
    sin_alpha = np.sin(np.deg2rad(alpha))

    a = phi * (phi / sigma * sin_alpha - J)
    b = sigma / 2 * J * P1 * (J * sin_alpha - phi)

    return a + b

def calc_C_thrust_y(phi, sigma, alpha, J, P1):
    cos_alpha = np.cos(np.deg2rad(alpha))

    a = phi**2 / sigma
    b = sigma / 2 * J**2 * P1

    return (a - b) * cos_alpha

def calc_C_power(phi, psi):
    return phi * psi

def calc_psi(phi, Cpr, J, eta):
    a = phi**2 * (1 - Cpr) - J**2

    return a / (2 * eta)