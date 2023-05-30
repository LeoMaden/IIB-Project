import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

plt.rcParams.update({"font.size": 20})

data = scipy.io.loadmat("CFD Data/diffuser_Cpr_xy_data.mat")["data"]

y = data[0, :, :]
z = data[1, :, :]
Cpr = data[2, :, :]

max_r = np.max(np.sqrt(y**2 + z**2))

fig, ax = plt.subplots(figsize=(10, 8))
CS = ax.contourf(z / max_r, y / max_r, Cpr, cmap="gray")
CS2 = ax.contour(CS, levels=[0], colors="r")
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel("$ C_{pr} $", rotation=0)
cbar.add_lines(CS2)
ax.set_aspect("equal")
ax.text(0.036 / max_r, 0.03 / max_r, "Suction surface", rotation=38, verticalalignment="top", rotation_mode="anchor")
ax.set_xlabel(r"$z/r_{max}$")
ax.set_ylabel(r"$\dfrac{y}{r_{max}}$", rotation=0, labelpad=20)
ax.text(0.05 / max_r, 0.018 / max_r, "$ \overline{C_{pr}} = " + str(-0.1895) + "$", bbox={"edgecolor": "k", "fill": None})

fig.savefig("Figures/Passage_diffuser_Cpr_contour.pdf")
plt.show()

