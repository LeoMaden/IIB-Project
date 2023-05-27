import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

from modules import experimental

plt.rcParams.update({"font.size": 20})


data = experimental.load_data(["Experimental Data/run3_processed.mat"])
markers = ["x", "o", "s", "+", "D", "^", "*"]

for run in data:

    R = 287
    rho = run.atm.Pa / (R * run.atm.Ta)

    V = np.sqrt(2 * (run.P_total - run.P_static) / rho)
    U = 2 * np.pi * run.f_shaft * run.geom.rmean
    rpm = 60 * run.f_shaft

    J = V / U

    A3 = run.geom.A1

    C_Fx = run.F_long / (rho * A3 * U**2)
    C_Fy = run.F_lat / (rho * A3 * U**2)

    alpha = run.aoa
    
    # Figure 1
    fig1, axs1 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
    labels1 = []
    handles1 = []

    # Figure 2
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    J_zero_C_Tx = np.zeros(len(alpha[:-1]))
    C_Fy_fit_params = np.zeros((len(alpha[:-1]), 3))

    # Figure 3
    fig3, ax3 = plt.subplots(figsize=(10, 8))


    # fig3, axs3 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
    # fig4, ax4 = plt.subplots(figsize=(10, 8))
    # fig5, axs5 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")

    # Fitting function
    fit_func = lambda x, a, b, c: a + b*x + c*x**2





    # --- Figure 1: C_Fx and C_Fy with quadratic fit lines vs J ---
    for i, alpha_i in enumerate(alpha[:-1]):

        cos_alpha_i = np.cos(np.deg2rad(alpha_i))
        sin_alpha_i = np.sin(np.deg2rad(alpha_i))

        # Pick out data at alpha_i
        # - Omit U=0 data
        # - Omit V=0 data
        J_i = J[i, 1:, 1:]
        C_Fx_i = C_Fx[i, 1:, 1:]
        C_Fy_i = C_Fy[i, 1:, 1:]
        C_F_i = C_Fy_i * cos_alpha_i + C_Fx_i * sin_alpha_i

        J_i_flat = J_i.flatten()
        C_Fx_i_flat = C_Fx_i.flatten()
        C_Fy_i_flat = C_Fy_i.flatten()
        C_F_i_flat = C_F_i.flatten()

        # Indices without nans
        i_C_Fx_valid = np.isfinite(J_i_flat) & np.isfinite(C_Fx_i_flat)
        i_C_Fy_valid = np.isfinite(J_i_flat) & np.isfinite(C_Fy_i_flat)
        i_C_F_valid = np.isfinite(J_i_flat) & np.isfinite(C_F_i_flat)

        # Fit quadratics
        min_J = 0
        max_J = max(J_i_flat[i_C_F_valid])

        J_i_linspace = np.linspace(min_J, max_J, 100)

        C_Fx_i_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i_flat[i_C_Fx_valid], C_Fx_i_flat[i_C_Fx_valid])
        C_Fy_i_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i_flat[i_C_Fy_valid], C_Fy_i_flat[i_C_Fy_valid])
        C_F_i_fit_params, _ = scipy.optimize.curve_fit(fit_func, J_i_flat[i_C_F_valid], C_F_i_flat[i_C_F_valid])

        C_Fx_i_fit = fit_func(J_i_linspace, *C_Fx_i_fit_params)
        C_Fy_i_fit = fit_func(J_i_linspace, *C_Fy_i_fit_params)
        C_F_i_fit = fit_func(J_i_linspace, *C_F_i_fit_params)

        color = list(mcolors.TABLEAU_COLORS.keys())[i]

        axs1[0].set_title("$ C_{F_x} $")
        axs1[0].plot(J_i_flat, C_Fx_i_flat, 'x', color=color, label=f"$ {alpha_i}^\circ $")
        axs1[0].plot(J_i_linspace, C_Fx_i_fit, '-', color=color)

        axs1[1].set_title("$ C_{F_y} $")
        axs1[1].plot(J_i_flat, C_Fy_i_flat, 'x', color=color, label=f"$ {alpha_i}^\circ $")
        axs1[1].plot(J_i_linspace, C_Fy_i_fit, '-', color=color)

        labels1.append(f"$ {alpha_i}^\circ $")
        handles1.append(mlines.Line2D([], [], linestyle='-', marker='x', color=color))

        C_Fy_fit_params[i] = C_Fy_i_fit_params

        # J where C_Tx = 0
        f = lambda J: fit_func(J, *C_Fx_i_fit_params)
        sol = scipy.optimize.root(f, 0)

        if sol.success == False:
            J_zero_C_Tx[alpha_i] = np.nan
            continue
            
        # Solution found
        J_zero_C_Tx[i] = sol.x[0]


    # --- Figure 2: C_Fx = 0 J vs. alpha ---

    # Logical to have J = 0 at alpha = 0
    J_zero_C_Tx[0] = 0

    ax2.plot(alpha[:-1], J_zero_C_Tx, 'x-')


    # --- Figure 3: 
    C_Fy_zero_C_Tx = np.zeros(len(alpha[:-1]))

    for i, alpha_i in enumerate(alpha[:-1]):
        
        # Read C_Ty at J where C_Tx = 0
        J_i_zero_C_Tx = J_zero_C_Tx[i]

        fit_params_i = C_Fy_fit_params[i]

        C_Fy_zero_C_Tx[i] = fit_func(J_i_zero_C_Tx, *fit_params_i)

    ax3.plot(alpha[:-1], C_Fy_zero_C_Tx, 'x-')



        # # --- Fig 2: Experiment C_F with quadratic fit lines vs J ---
        # ax2.plot(J_i_flat, C_F_i_flat, 'x', color=color, label=alpha_i)
        # ax2.plot(J_i_linspace, C_F_i_fit, '-', color=color)

        # # --- Model evaluation ---
        # phi = 0.6
        # sigma = 1.3
        # Mb_model = 0.05
        # gam = 1.4

        # C_Tx_model = phi * (phi / sigma * sin_alpha_i - J_i_linspace)
        # C_Ty_model = np.ones_like(J_i_linspace) * phi**2 / sigma * cos_alpha_i
        # C_T_model = C_Ty_model * cos_alpha_i + C_Tx_model * sin_alpha_i

        # # --- Fig 3: Model C_Tx and C_Ty and experiment C_Fx and C_Fy vs J ---
        # axs3[0].plot(J_i_linspace, C_Tx_model, '-', color=color, label=alpha_i)
        # axs3[0].plot(J_i_linspace, C_Fx_i_fit, '--', color=color)

        # axs3[1].plot(J_i_linspace, C_Ty_model, '-', color=color, label=alpha_i)
        # axs3[1].plot(J_i_linspace, C_Fy_i_fit, '--', color=color)

        # # --- Fig 4: Model C_T and experiment C_F vs J ---
        # ax4.plot(J_i_linspace, C_T_model, '-', color=color, label=alpha_i)
        # ax4.plot(J_i_linspace, C_F_i_fit, '--', color=color)

        # # --- Fig 5: Difference between model C_Ti and experiment C_Fi ---
        # C_Fx_i_trunc_fit = C_Fx_i_fit_func.convert().cutdeg(1)(J_i_linspace)
        # C_Fy_i_trunc_fit = C_Fy_i_fit_func.convert().cutdeg(1)(J_i_linspace)

        # # axs5[0].plot(J_i_linspace, C_Fx_i_fit, '--', color=color, label=alpha_i)
        # axs5[0].plot(J_i_linspace, C_Fx_i_trunc_fit, '--', color=color, label=alpha_i)
        # axs5[0].plot(J_i_linspace, C_Tx_model, '-', color=color, label=alpha_i)

        # # axs5[1].plot(J_i_linspace, C_Fy_i_fit, '--', color=color, label=alpha_i)
        # axs5[1].plot(J_i_linspace, C_Fy_i_trunc_fit, '--', color=color, label=alpha_i)
        # axs5[1].plot(J_i_linspace, C_Ty_model, '-', color=color, label=alpha_i)

        # # axs5[0].plot(J_i_linspace, C_Fx_i_fit - C_Tx_model, '-', color=color, label=alpha_i)
        # # axs5[1].plot(J_i_linspace, C_Fy_i_fit - C_Ty_model, '-', color=color, label=alpha_i)


# Figure 1 properties
for ax in axs1:
    ax.set_xlabel("$ J $")
    ax.grid('major')

axs1[1].legend(labels=labels1, handles=handles1, title=r"$ \alpha $")

fig1.savefig("Figures/Exp_LiftFan_xy_forces.pdf")

# Figure 2 properties
ax2.set_title(r"$ C_{F_x} = 0 $")
ax2.set_xlabel(r"$ \alpha $ (deg)")
ax2.set_ylabel(r"$ J $ ", rotation=0, labelpad=20)
ax2.grid('major')

fig2.savefig("Figures/Exp_LiftFan_alpha_J_zero_x_force.pdf")


# Figure 3 properties
ax3.set_title(r"$ C_{F_x} = 0 $")
ax3.set_xlabel(r"$ \alpha $ (deg)")
ax3.set_ylabel(r"$ C_{F_y} $", rotation=0, labelpad=20)
ax3.grid('major')

fig3.savefig("Figures/Exp_LiftFan_alpha_C_Ty_zero_x_force.pdf")


# ax2.legend()
# axs3[0].legend()
# axs3[1].legend()
# ax4.legend()
plt.show()
