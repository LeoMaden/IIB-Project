import matplotlib.pyplot as plt
import numpy as np
from modules import aero_model


sigma = np.linspace(1, 2, 100)
Cpr = 1 - 1/sigma**2

eta = 0.9
J = 0.25
phi = 0.7

psi = aero_model.calc_psi(phi, Cpr, J, eta)
C_Wdotx = aero_model.calc_C_power(phi, psi)

plt.plot(sigma, J/C_Wdotx)
plt.show()