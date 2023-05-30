import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker

plt.rcParams.update({"font.size": 20})


# Load rotor flow data
col_names = ["r_guess", "alpha_guess", "To_guess", "Po_guess", "r_CFD", "alpha_CFD", "To_CFD", "Po_CFD"]
data = pd.read_csv("CFD Data/rotor_outlet_flow.csv", names=col_names)

# # Convert to kPa
# data["Po_guess"] /= 1000
# data["Po_real"] /= 1000

# Plot data
# fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(15, 8), sharey="row")
fig, ax = plt.subplots(figsize=(10, 8))

ax.plot("alpha_guess", "r_guess", data=data)
ax.plot("alpha_CFD", "r_CFD", data=data)
ax.set_xlabel(r"$ \alpha $ (deg)")
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(15))
ax.set_xlim(xmin=0)
ax.set_ylabel(r"$ \dfrac{r - r_{hub}}{r_{cas} - r_{hub}} $", rotation=90)
ax.set_ylim((0, 1))
ax.legend(["Guess", "CFD"])

# axs[1].plot("To_guess", "r_guess", data=data)
# axs[1].plot("To_real", "r_real", data=data)
# axs[1].set_xlabel(r"$ T_0 $ (K)")

# axs[2].plot("Po_guess", "r_guess", data=data, label="Guess")
# axs[2].plot("Po_real", "r_real", data=data, label="Real")
# axs[2].set_xlabel(r"$ P_0 $ (kPa)")

# axs[2].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#            ncol=2, mode="expand", borderaxespad=0.)
# axs[2].legend(bbox_to_anchor=(0, 0.92), loc="upper left")

# for ax in axs:
ax.grid("major")

fig.set_constrained_layout(True)
fig.savefig("Figures/Rotor_outlet_flow.pdf")
plt.show()