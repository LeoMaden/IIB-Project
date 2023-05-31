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
fig, ax = plt.subplots(figsize=(10, 8))
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
    # labels += [fr"{rpm_i}, ${alpha_j}^\circ $" for alpha_j in alpha]
    labels += [fr"${alpha_j}^\circ $" for alpha_j in alpha]

    # Calculate pressure coefficient
    Cp = (Pe - P1) / (0.5 * rho * V**2)

    # Calculate Reynolds number
    ReD = rho * V * D / mu

    # Blade Reynolds number
    U = 2 * np.pi * Rmean * rpm_i / 60
    ReU = np.mean(rho) * U * D / mu

    J = V / U

    for j, alpha_j in enumerate(alpha):
        V_j = V[:, j]
        J_j = J[:, j]
        ReD_j = ReD[:, j]
        Cp_j = Cp[:, j]

        # ax.loglog(ReD_j, Cp_j, marker=markers[i], linestyle='', color=colors[j])
        ax.plot(1/J_j, Cp_j, marker=markers[i], linestyle='', color=colors[j])

    # Trendline
    fit_func = lambda x, a, b: a/(x**b)

    # ReD_flat = ReD.flatten()
    J_flat = J.flatten()

    Cp_flat = Cp.flatten()

    # i_valid = np.isfinite(ReD_flat) & np.isfinite(Cp_flat) 
    i_valid = np.isfinite(J_flat) & np.isfinite(Cp_flat) 

    # ReD_fit_params, _ = scipy.optimize.curve_fit(fit_func, ReD_flat[i_valid], Cp_flat[i_valid])
    J_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_flat[i_valid], Cp_flat[i_valid])

    # ReD_fit = np.linspace(min(ReD_flat[i_valid]), (max(ReD_flat[i_valid])), 100)
    J_fit = np.linspace(min(J_flat[i_valid]), (max(J_flat[i_valid])), 100)

    # Cp_fit = fit_func(ReD_fit, *ReD_fit_params)
    Cp_fit = fit_func(J_fit, *J_fit_params)

    # ax.loglog(ReD_fit, Cp_fit, 'k-')
    ax.plot(1/J_fit, Cp_fit, 'k-')


# ax.set_ylim(-30, 30)
ax.grid('major')
# ax.legend(labels, title=r"rpm, $ \alpha $", bbox_to_anchor=(1.05, 1), loc='upper left')
ax.legend(labels, title=r"$ \alpha $")
ax.set_ylabel(r"$ \dfrac{P_e - P_1}{\frac{1}{2} \rho V^2} $", rotation=0, labelpad=40)

# ax.set_xlabel(r"$Re_D$")
ax.set_xlabel(r"$1/J$")

ax.set_title(f"$ Re_U = {ReU:.3g} $")

fig.set_constrained_layout(True)
fig.savefig("Figures/IPM5_CFD_exit_pressure.pdf")
plt.show()
