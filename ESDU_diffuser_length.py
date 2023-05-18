"""
Look at how the diffuser length for optimum pressure recovery
changes as the area ratio is increased 
"""
from modules import esdu_75026 as esdu
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})

# Data for area ratio and optimum pressure recovery
# from EDSU
A2_A1 = esdu.Cpr2_data_symmetrical["A2/A1"].to_numpy()
L_h1 = esdu.Cpr2_data_symmetrical["L/h1"].to_numpy()

poly_fit = np.polynomial.polynomial.Polynomial.fit(A2_A1, L_h1, deg=2)
fit_x = np.linspace(1.1, 3.5, 100)
fit_y = poly_fit(fit_x)
a, b, c = poly_fit.convert().coef
fit_str = fr"$ {a:.3g} + {b:.3g}x + {c:.3g}x^2 $"

# Plot results
plt.figure(figsize=(10, 8))
plt.plot(A2_A1, L_h1, marker="x", label="Data")
plt.plot(fit_x, fit_y, label=f"Fit: {fit_str}")

plt.xlabel(r"$ \dfrac{A_2}{A_1} $")
plt.ylabel(r"$ \dfrac{L}{h_1} $", rotation=0, labelpad=20)
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig("Figures/ESDU_diffuser_length.pdf")
plt.show()