"""
Compare model and experimental data when the wind tunnel is
switched off (J=0) and the load cell is reading the thrust
components directly
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from modules import aero_model
from modules import experimental

plt.rcParams.update({"font.size": 20})

# LiftFan parameters
sigma = 1.3
phi = 0.6
J = 0
P1 = 0

# Model values
alpha = np.linspace(0, 90, 100)
C_Tx = aero_model.calc_C_thrust_x(phi, sigma, alpha, J, P1)
C_Ty = aero_model.calc_C_thrust_y(phi, sigma, alpha, J, P1)

# --- Experimental values ---
# Load data
runs = []
base_path = "./Experimental Data/{}_processed.mat"
paths = [base_path.format(run) for run in runs]
data = experimental.load_data(paths)

# TODO: Analyse data


# Plot
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(12, 8))

axs[0].plot(alpha, C_Tx)
axs[0].set_title("$ C_{T_x} $")

axs[1].plot(alpha, C_Ty)
axs[1].set_title("$ C_{T_y} $")

for ax in axs:
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(45))
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(15))
    ax.set_xlabel(r"$\alpha$ (deg)")
    ax.grid(which="both")

plt.show()