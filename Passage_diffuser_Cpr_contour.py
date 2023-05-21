import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

plt.rcParams.update({"font.size": 20})

data = scipy.io.loadmat("CFD Data/diffuser_Cpr_xy_data.mat")["data"]

x = data[0, :, :]
y = data[1, :, :]
Cpr = data[2, :, :]

fig, ax = plt.subplots(figsize=(10, 8))
CS = ax.contourf(x, y, Cpr, cmap="gray")
CS2 = ax.contour(CS, levels=[0], colors="r")
cbar = fig.colorbar(CS)
cbar.ax.set_ylabel("$ C_{pr} $", rotation=0)
cbar.add_lines(CS2)
ax.set_aspect("equal")
ax.text(0.03, 0.036, "Suction surface", rotation=52, verticalalignment="bottom", rotation_mode="anchor")
ax.set(xlabel="$x$", ylabel="$y$")
ax.text(0.02, 0.08, "$ \overline{C_{pr}} = " + str(-0.1895) + "$", bbox={"color": "grey"})

fig.savefig("Figures/Passage_diffuser_Cpr_contour.pdf")
plt.show()

