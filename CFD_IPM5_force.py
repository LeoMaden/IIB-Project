import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy.optimize

plt.rcParams.update({"font.size": 20})

# Parameters
rpm = 5000
D = 0.16
Rc = D / 2
hub_tip = 0.35
Rh = hub_tip * Rc
R_mean = 0.5 * (Rc + Rh)
U = 2 * np.pi * R_mean * rpm / 60
A3 = np.pi * (Rc**2 - Rh**2)

print(f"U = {U}")

# Load data
filename = f"CFD Data/IPM5_force_{rpm}rpm.mat"
data = scipy.io.loadmat(filename)

V = data["V"]
rho = data["rho"]
Fx = data["Fx"]
Fy = data["Fy"]
Vf = data["Vf"][0, :]
alpha = data["alpha"][0, :]
mask = data["mask"] == 1

# Calculate nondimensionals
J = V / U
C_Fx = Fx / (rho * A3 * U**2)
C_Fy = Fy / (rho * A3 * U**2)

# Plotting 
fig1, axs1 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
colors = list(mcolors.TABLEAU_COLORS.keys())

# Fitting
fit_func = lambda x, a, b, c: a*x**2 + b*x + c

for i, alpha_i in enumerate(alpha):
    J_i = J[:, i]
    C_Fx_i = C_Fx[:, i]
    C_Fy_i = C_Fy[:, i]

    i_valid = np.isfinite(J_i) & np.isfinite(C_Fx_i) & np.isfinite(C_Fy_i) & mask[:, i]

    J_i = J_i[i_valid]
    C_Fx_i = C_Fx_i[i_valid]
    C_Fy_i = C_Fy_i[i_valid]

    J_linspace = np.linspace(0, max(J_i), 100)

    C_Fx_i_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i, C_Fx_i)
    C_Fy_i_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i, C_Fy_i)

    C_Fx_i_fit = fit_func(J_linspace, *C_Fx_i_fit_params)
    C_Fy_i_fit = fit_func(J_linspace, *C_Fy_i_fit_params)

    axs1[0].plot(J_i, C_Fx_i, 'x', color=colors[i])
    axs1[0].plot(J_linspace, C_Fx_i_fit, '-', color=colors[i])

    axs1[1].plot(J_i, C_Fy_i, 'x', color=colors[i])
    axs1[1].plot(J_linspace, C_Fy_i_fit, '-', color=colors[i])

# Figure 1 settings
axs1[0].set_title(r"$ C_{T_x} $")
axs1[1].set_title(r"$ C_{T_y} $")

for ax in axs1:
    ax.set_xlabel(r"$ J $")
    ax.grid('major')

fig1.set_constrained_layout(True)

plt.show()



