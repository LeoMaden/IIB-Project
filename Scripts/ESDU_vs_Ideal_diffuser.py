"""
Look at the pressure recovery for a symmetrical annular diffuser
using data from ESDU 75026 and compare to an ideal diffuser
"""
import sys
sys.path.append("./Modules")

import esdu_75026 as esdu
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

plt.rcParams.update({"font.size": 20})

# Data for area ratio and optimum pressure recovery
# from EDSU
A2_A1 = esdu.Cpr2_data_symmetrical["A2/A1"]
Cpr_opt = esdu.Cpr2_data_symmetrical["Cpr"]

# Ideal diffuser using Bernoulli
A2_A1_ideal = np.linspace(min(A2_A1), max(A2_A1), 100)
Cpr_opt_ideal = 1 - 1 / A2_A1_ideal**2

# ESDU Fit
fit_func = lambda x, a, b: a - b / x**2
(a_fit, b_fit), _ = scipy.optimize.curve_fit(fit_func, A2_A1, Cpr_opt)
fit_x = A2_A1_ideal
fit_y = fit_func(fit_x, a_fit, b_fit)
fit_str = r"$" + f"{a_fit:.3g} - {b_fit:.3g}/x^2 $"# + r"(A_2/A_1)^{-2} $"

# Plot results
plt.figure(figsize=(10, 8))
plt.plot(A2_A1_ideal, Cpr_opt_ideal, label="Ideal")
plt.plot(A2_A1, Cpr_opt, marker="x", label="ESDU")
plt.plot(fit_x, fit_y, label=f"Fit: {fit_str}")

plt.xlabel(r"$ \dfrac{A_2}{A_1} $")
plt.ylabel(r"$ C_{pr}^{**} $", rotation=0, labelpad=20)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("Figures/ESDU_vs_Ideal_diffuser.pdf")
plt.show()



