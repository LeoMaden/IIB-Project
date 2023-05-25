import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import matplotlib.colors as mcolors

plt.rcParams.update({"font.size": 20})

rpm = [22000]

# Air viscosity
mu = 1.81e-5

# Geometry
D = 0.16

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
    mask = data["mask"]

    i_1 = list(stations).index("Domain inlet")
    i_e = list(stations).index("Diffuser inlet")

    V = data["V"][:, :, i_1]
    P1 = data["P"][:, :, i_1]
    rho = data["ro"][:, :, i_1]
    Pe = data["P"][:, :, i_e]

    # Plotting setup
    labels += [fr"{rpm_i}, ${alpha_j}^\circ $" for alpha_j in alpha]

    # Calculate pressure coefficient
    Cp = (Pe - P1) / (0.5 * rho * V**2)

    # Calculate Reynolds number
    ReD = rho * V * D / mu

    for j, alpha_j in enumerate(alpha):
        V_j = V[:, j]
        ReD_j = ReD[:, j]
        Cp_j = Cp[:, j]
        ax.loglog(ReD_j, Cp_j, marker=markers[i], linestyle='', color=colors[j])

    # Trendline
    fit_func = lambda x, a, b: a/(x**b)

    ReD_flat = ReD.flatten()
    Cp_flat = Cp.flatten()  
    i_valid = np.isfinite(ReD_flat) & np.isfinite(Cp_flat)

    fit_params, _ = scipy.optimize.curve_fit(fit_func, ReD_flat[i_valid], Cp_flat[i_valid])
    print(fit_params)

    ReD_fit = np.linspace(min(ReD_flat[i_valid]), (max(ReD_flat[i_valid])), 100)
    Cp_fit = fit_func(ReD_fit, *fit_params)

    ax.loglog(ReD_fit, Cp_fit, 'k-')


# ax.set_ylim(-30, 30)
ax.grid('major')
ax.legend(labels, title=r"rpm, $ \alpha $", bbox_to_anchor=(1.05, 1), loc='upper left')
ax.set_ylabel(r"$ \dfrac{P_e - P_1}{\frac{1}{2} \rho V^2} $", rotation=0, labelpad=30)
ax.set_xlabel(r"$Re_D$")
fig.set_constrained_layout(True)
fig.savefig("Figures/IPM5_CFD_exit_pressure.pdf")
plt.show()
