"""
Look at how the diffuser length for optimum pressure recovery
changes as the area ratio is increased 
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
L_h1 = esdu.Cpr2_data_symmetrical["L/h1"]

# Plot results
plt.figure(figsize=(10, 8))
plt.plot(sigma, L_h1, marker="x")

plt.xlabel(r"$ \sigma $")
plt.ylabel(r"$ \dfrac{L}{h_1} $", rotation=0, labelpad=20)
plt.grid()
plt.tight_layout()
plt.savefig("Figures/ESDU_diffuser_length.pdf")
plt.show()