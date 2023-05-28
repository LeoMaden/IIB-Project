import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.size": 20})


# Geometry
hub_tip = 0.3
D2_A3 = 1 / (4 * np.pi) / (1 - hub_tip**2)

# Drag
f = 1.28

# --- Figure 1: C_Dx vs alpha and L/D ---
L_D = np.linspace(0.3, 3, 100)
alpha = np.linspace(0, 50, 100)

_L_D, _alpha = np.meshgrid(L_D, alpha)

_sin_alpha = np.sin(np.deg2rad(_alpha))
_cos_alpha = np.cos(np.deg2rad(_alpha))

_C_Dx = f * D2_A3 * (0.25 * np.pi * _sin_alpha + _L_D * _cos_alpha)

# Plot
fig1, ax1 = plt.subplots(figsize=(10, 8))
con = ax1.contourf(_alpha, _C_Dx, _L_D, cmap="viridis", levels=100)

cbar = fig1.colorbar(con, ax=ax1)
cbar.set_label(r"$ \dfrac{L}{D} $", rotation=0, labelpad=20)
ax1.set_xlabel(r"$ \alpha $ (deg)")
ax1.set_ylabel(r"$ C_{D_x} $", rotation=0, labelpad=20)

fig1.savefig("Figures/Model_drag_x.pdf")

plt.show()



