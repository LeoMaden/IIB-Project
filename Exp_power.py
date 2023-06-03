import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from modules import aero_model

from modules import experimental

plt.rcParams.update({"font.size": 20})


data = experimental.load_data(["Experimental Data/run3_processed.mat"])[0]


R = 287
rho = data.atm.Pa / (R * data.atm.Ta)

V = np.sqrt(2 * (data.P_total - data.P_static) / rho)
U = 2 * np.pi * data.f_shaft * data.geom.rmean
rpm = 60 * data.f_shaft
Wdotx = data.V_supply * data.I_supply
A3 = data.geom.A1
alpha = data.aoa

J = V / U
C_Wdotx = Wdotx / (rho * A3 * U**3)


colors = list(mcolors.TABLEAU_COLORS.keys())

# Figure 1
fig1, ax1 = plt.subplots(figsize=(10, 8))
labels1 = []
handles1 = []

fig2, ax2 = plt.subplots(figsize=(10, 8))
labels2 = []
handles2 = []

fit_func = lambda x, a, b, c: a + b*x + c*x**2

for i, alpha_i in enumerate(alpha):
    if alpha_i not in [10, 40, 90]: continue
    # if alpha_i == 90: continue

    # Data for pwm speed 2: as current zero before that
    # Ignore zero tunnel speed
    J_i = J[i, :, 2:].flatten()
    C_Wdotx_i = C_Wdotx[i, :, 2:].flatten()

    i_valid = np.isfinite(J_i) & np.isfinite(C_Wdotx_i) 

    J_i = J_i[i_valid]
    C_Wdotx_i = C_Wdotx_i[i_valid]

    fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i, C_Wdotx_i)
    _J = np.linspace(0, max(J_i))
    _C_Wdotx = fit_func(_J, *fit_params)

    ax1.plot(J_i, C_Wdotx_i, 'x', color=colors[i])
    ax1.plot(_J, _C_Wdotx, '-', color=colors[i])

    handles1.append(mlines.Line2D([], [], linestyle='-', marker='x', color=colors[i]))
    labels1.append(fr"$ {alpha_i}^\circ $")


    if alpha_i in [10, 40]:
        h, = ax2.plot(_J, _C_Wdotx, ':', color=colors[i])
        handles2.append(h)    
        labels2.append(fr"$ \alpha = {alpha_i}^\circ $")


# Plot model
phi = 0.6
sigma = 1.3
eta = 0.9
Cpr = 1 - 1/sigma**2
# Cpr = -1
_J = np.linspace(0, 0.7, 100)
C_Wdotx_model = phi / (2 * eta) * (phi**2 * (1 - Cpr) - _J**2)

ax1.plot(_J, C_Wdotx_model, '-k')
model_line = mlines.Line2D([], [], linestyle='-', color='k')
fig1.add_artist(ax1.legend([model_line], ["Model"], loc='lower left'))

ax1.set_ylim(bottom=0)
ax1.legend(handles1, labels1, title=r"$ \alpha $", loc='lower right')
ax1.set_xlabel(r"$ J $")
ax1.set_ylabel(r"$ C_{\dot W_x} $", rotation=0, labelpad=20)
ax1.grid('major')
fig1.set_constrained_layout(True)
fig1.savefig("Figures/Exp_power.pdf")


# Add legend for experimental data
# fig2.add_artist(ax2.legend(handles2, labels1, title=r"$ \alpha $", loc='upper right', ncols=2))

# Old and new model plot
new_colors = ['fuchsia', 'lime']
new_colors = ['tab:purple', 'tab:orange']
Ypi = [1.1, 1.2]
# labels2 = []
# handles2 = []

for color, Ypi_i in zip(new_colors, Ypi):
    C_Wdotx_new_model1 = phi / (2 * eta) * (phi**2 * (1 - Cpr) - _J**2 * (1-Ypi_i))
    h, = ax2.plot(_J, C_Wdotx_new_model1, color=color, linewidth=2)
    handles2.append(h)
    labels2.append("$ Y_{p,i} = " + str(Ypi_i) + "$")

fig2.add_artist(ax2.legend(handles2, labels2))


# ax2.set_ylim(bottom=0)
ax2.set_xlabel(r"$ J $")
ax2.set_ylabel(r"$ C_{\dot W_x} $", rotation=0, labelpad=20)
ax2.grid('major')
fig2.set_constrained_layout(True)

fig2.savefig("Figures/Exp_power_new_model.pdf")


plt.show()
