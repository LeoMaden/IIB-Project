"""
Look at the pressure recovery for a symmetrical annular diffuser
using data from ESDU 75026 and compare to an ideal diffuser
"""
import sys
sys.path.append("./Modules")

import esdu_75026 as esdu
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})

# Data for area ratio and optimum pressure recovery
# from EDSU
sigma = esdu.Cpr2_data_symmetrical["A2/A1"]
Cpr_opt = esdu.Cpr2_data_symmetrical["Cpr"]

# Ideal diffuser using Bernoulli
sigma_ideal = np.linspace(min(sigma), max(sigma), 100)
Cpr_opt_ideal = 1 - 1 / sigma_ideal**2

# Plot results
plt.figure(figsize=(10, 8))
plt.plot(sigma_ideal, Cpr_opt_ideal, label="Ideal")
plt.plot(sigma, Cpr_opt, marker="x", label="ESDU")

plt.xlabel(r"$ \dfrac{A_2}{A_1} $")
plt.ylabel(r"$ C_{pr}^{**} $")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Figures/ESDU_vs_Ideal_diffuser.pdf")
plt.show()



