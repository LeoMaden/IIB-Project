import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import matplotlib.colors as mcolors

plt.rcParams.update({"font.size": 20})

rpm = [22000]

# Air viscosity
mu = 1.81e-5

# Design params
D = 0.16
hub_tip = 0.35
phi_des = 0.7

Rc = D / 2
Rh = Rc * hub_tip
Rmean = 0.5 * (Rc + Rh)

# Setup plot
fig, ax = plt.subplots(figsize=(10, 7))
labels = []
markers = ['x', 'o']
colors = list(mcolors.TABLEAU_COLORS.keys())

for i, rpm_i in enumerate(rpm):

    # Load data
    filename = f"CFD Data/IPM5_avg_flow_properties_{rpm_i}rpm.mat"
    data = scipy.io.loadmat(filename, squeeze_me=True, chars_as_strings=True)

    # Extract variables
    Vf = data["Vf"]
    alpha = data["alpha"]
    stations = data["stations"]
    mask = data["mask"] == 1

    i_1 = list(stations).index("Domain inlet")
    i_3 = list(stations).index("Rotor inlet")

    V = data["V"][:, :, i_1]
    V3 = data["V"][:, :, i_3]
    U = 2 * np.pi * Rmean * rpm_i / 60

    # Plotting setup
    labels += [fr"{rpm_i}, ${alpha_j}^\circ $" for alpha_j in alpha]

    # Calculate flow coefficient
    phi = V3/U

    # Calculate Advance ratio
    J = V / U

    # Add nans at invalid values
    J[~mask] = np.nan

    for j, alpha_j in enumerate(alpha):


        phi_j = phi[:, j]
        J_j = J[:, j]
        ax.plot(J_j, phi_j, marker=markers[i], linestyle='', color=colors[j])
        

# Plot design phi
phi_des_line = ax.hlines(phi_des, min(J.flatten()), max(J.flatten()), color='k', linestyle='--')
legend = ax.legend(handles=[phi_des_line], labels=["Design $\phi$"], loc='lower right')
ax.add_artist(legend)
ax.set_ylim(0, 1.1*phi_des)
ax.set_xlim(left=0)



ax.grid('major')
ax.legend(labels, title=r"rpm, $ \alpha $", bbox_to_anchor=(1, 1), loc='upper left')
ax.set_ylabel(r"$ \phi $", rotation=0, labelpad=15)
ax.set_xlabel(r"$J$")
fig.set_constrained_layout(True)
fig.savefig("Figures/IPM5_CFD_flow_coefficient.pdf")
plt.show()
