from modules import experimental
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import scipy.optimize
import scipy.interpolate
import matplotlib.lines as mlines

plt.rcParams.update({"font.size": 20})

# Load data
data = experimental.load_data(["Experimental Data/run3_processed.mat"])

# Define figures
fig1, axs1 = plt.subplots(nrows=1, ncols=2, figsize=(12, 8), sharex='all', sharey='all')
fig2, axs2 = plt.subplots(nrows=1, ncols=2, figsize=(12, 8), sharex='all', sharey='all')

# Look at case where fan is off (U=0) and tunnel is on (V>0)
for run in data:
    i_fan_off = list(run.pwm_speed).index(0)
    i_tunnel_on = (run.tunnel_switch != 0)

    P_static = run.P_static[:, i_tunnel_on, i_fan_off]
    P_total = run.P_total[:, i_tunnel_on, i_fan_off]

    F_x = run.F_long[:, i_tunnel_on, i_fan_off]
    F_y = run.F_lat[:, i_tunnel_on, i_fan_off]

    R = 287
    mu = 1.81e-5
    rho = run.atm.Pa / (R * run.atm.Ta)

    V = np.sqrt(2 * (P_total - P_static) / rho)

    A3 = run.geom.A1
    D = 2 * run.geom.rtip

    ReD = rho * V * D / mu

    C_Dx = -F_x / (0.5 * rho * A3 * V**2)
    C_Dy = -F_y / (0.5 * rho * A3 * V**2)

    alpha = run.aoa

    colors = list(mcolors.TABLEAU_COLORS.keys())
    leg_handles1 = []
    leg_labels1 = []

    for i, alpha_i in enumerate(alpha):
        V_i = V[i, :]
        ReD_i = ReD[i, :]
        C_Dx_i = C_Dx[i, :]
        C_Dy_i = C_Dy[i, :]

        # ReD_i_fit = np.linspace(min(ReD_i), max(ReD_i), 100)
        ReD_i_fit = np.linspace(0.5e5, 4e5, 100)

        xfit_func = lambda x, a, b: a / x + b / np.sqrt(x)
        yfit_func = lambda x, a, b: a / x + b / np.sqrt(x)


        # C_Dx_i_fit_func = np.polynomial.Polynomial.fit(ReD_i, C_Dx_i, deg=2)
        # C_Dy_i_fit_func = np.polynomial.Polynomial.fit(ReD_i, C_Dy_i, deg=2)

        C_Dx_i_fit_params, _ = scipy.optimize.curve_fit(xfit_func, ReD_i, C_Dx_i)
        C_Dy_i_fit_params, _ = scipy.optimize.curve_fit(yfit_func, ReD_i, C_Dy_i)

        # C_Dx_i_fit = C_Dx_i_fit_func(ReD_i_fit)
        # C_Dy_i_fit = C_Dy_i_fit_func(ReD_i_fit)

        C_Dx_i_fit = xfit_func(ReD_i_fit, *C_Dx_i_fit_params)
        C_Dy_i_fit = yfit_func(ReD_i_fit, *C_Dy_i_fit_params)

        color = colors[i]

        axs1[0].semilogx(ReD_i, C_Dx_i, 'x', color=color)
        axs1[0].semilogx(ReD_i_fit, C_Dx_i_fit, '-', color=color, label=alpha_i)

        axs1[1].semilogx(ReD_i, C_Dy_i, 'x', color=color)
        axs1[1].semilogx(ReD_i_fit, C_Dy_i_fit, '-', color=color, label=alpha_i)

        leg_labels1.append(f"{alpha_i}°")
        leg_handles1.append(mlines.Line2D([], [], marker='x', linestyle='-', color=color))


    axs1[0].set_title(r"$ C_{D_x} $")
    axs1[1].set_title(r"$ C_{D_y} $")

    for ax in axs1:
        ax.tick_params('x', which='both', labelrotation=45)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(1e5))
        ax.xaxis.set_minor_locator(mticker.MultipleLocator(1e4))
        ax.xaxis.set_minor_formatter("")
        ax.grid('major')
    # axs1[1].xaxis.set_major_locator(mticker.MultipleLocator(1e5))
    # axs1[1].xaxis.set_major_formatter("")
    #axs1[1].xaxis.get_major_formatter().set_powerlimits((0, 3))
    #axs1[1].xaxis.get_major_formatter().set_scientific(True)
    #axs1[1].xaxis.get_major_formatter().set_useMathText(True)

    for ax in axs1:
        ax.set_xlabel("$ Re_D $")

    axs1[1].legend(labels=leg_labels1, handles=leg_handles1, bbox_to_anchor=(1.02, 1), loc='upper left', title=r"$ \alpha $")

    # --- Fig 2: Plot C_Dx and C_Dy vs alpha ---

    # C_D at highest ReD
    i_max_ReD = np.expand_dims(np.argmax(ReD, axis=1), axis=1)
    ReD_max = np.mean(np.take_along_axis(ReD, i_max_ReD, 1).squeeze())
    C_Dx_max_ReD = np.take_along_axis(C_Dx, i_max_ReD, 1).squeeze()
    C_Dy_max_ReD = np.take_along_axis(C_Dy, i_max_ReD, 1).squeeze()

    # C_D at lowest ReD
    i_min_ReD = np.expand_dims(np.argmin(ReD, axis=1), axis=1)
    ReD_min = np.mean(np.take_along_axis(ReD, i_min_ReD, 1).squeeze())
    C_Dx_min_ReD = np.take_along_axis(C_Dx, i_min_ReD, 1).squeeze()
    C_Dy_min_ReD = np.take_along_axis(C_Dy, i_min_ReD, 1).squeeze()

    # # Spline curves between points
    # alpha_linspace = np.linspace(0, 90, 91)

    # C_Dx_spline_func = scipy.interpolate.make_interp_spline(alpha, C_Dx)
    # C_Dy_spline_func = scipy.interpolate.make_interp_spline(alpha, C_Dy)

    # C_Dx_spline = C_Dx_spline_func(alpha_linspace)
    # C_Dy_spline = C_Dy_spline_func(alpha_linspace)

    leg_labels2 = [f"{ReD_max:.2e}", f"{ReD_min:.2e}"]

    axs2[0].set_title(r"$ C_{D_x} $")
    axs2[0].plot(alpha, C_Dx_max_ReD, 'x-', color='tab:blue')
    axs2[0].plot(alpha, C_Dx_min_ReD, 'x-', color='tab:orange')
    # axs2[0].plot(alpha_linspace, C_Dx_spline, '-', color='tab:blue')

    axs2[1].set_title(r"$ C_{D_y} $")
    axs2[1].plot(alpha, C_Dy_max_ReD, 'x-', color='tab:blue')
    axs2[1].plot(alpha, C_Dy_min_ReD, 'x-', color='tab:orange')
    # axs2[1].plot(alpha_linspace, C_Dy_spline, '-', color='tab:blue')

    for ax in axs2:
        ax.set_xlabel(r"$ \alpha $ (deg)")
        ax.grid('major')
    axs2[1].legend(title='Mean $ Re $', labels=leg_labels2)



fig1.set_constrained_layout(True)
fig2.tight_layout()

fig1.savefig("Figures/Exp_drag_coef_vs_ReD.pdf")
fig2.savefig("Figures/Exp_drag_coef_vs_alpha.pdf")

plt.show()