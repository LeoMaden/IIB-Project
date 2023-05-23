from modules import experimental
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

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
    
    # --- C_T vs J for specific angle ---
    fig1, axs1 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    fig3, axs3 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
    fig4, ax4 = plt.subplots(figsize=(10, 8))
    fig5, axs5 = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")

    for i, alpha_i in enumerate(alpha):
        cos_alpha_i = np.cos(np.deg2rad(alpha_i))
        sin_alpha_i = np.sin(np.deg2rad(alpha_i))

        J_i = J[i, :, 1:]
        C_Fx_i = C_Fx[i, :, 1:]
        C_Fy_i = C_Fy[i, :, 1:]
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
        min_J = min(J_i_flat[i_C_F_valid])
        max_J = max(J_i_flat[i_C_F_valid])

        J_i_linspace = np.linspace(min_J, max_J, 100)

        C_Fx_i_fit_func = np.polynomial.Polynomial.fit(J_i_flat[i_C_Fx_valid], C_Fx_i_flat[i_C_Fx_valid], deg=2)
        C_Fy_i_fit_func = np.polynomial.Polynomial.fit(J_i_flat[i_C_Fy_valid], C_Fy_i_flat[i_C_Fy_valid], deg=2)
        C_F_i_fit_func = np.polynomial.Polynomial.fit(J_i_flat[i_C_F_valid], C_F_i_flat[i_C_F_valid], deg=2)

        C_Fx_i_fit = C_Fx_i_fit_func(J_i_linspace)
        C_Fy_i_fit = C_Fy_i_fit_func(J_i_linspace)
        C_F_i_fit = C_F_i_fit_func(J_i_linspace)

        color = list(mcolors.TABLEAU_COLORS.keys())[i]

        # --- Fig 1: Experiment C_Fx and C_Fy with quadratic fit lines vs J ---
        axs1[0].set_title("$ C_{F_x} $")
        axs1[0].plot(J_i_flat, C_Fx_i_flat, 'x', color=color, label=alpha_i)
        axs1[0].plot(J_i_linspace, C_Fx_i_fit, '-', color=color)

        axs1[1].set_title("$ C_{F_y} $")
        axs1[1].plot(J_i_flat, C_Fy_i_flat, 'x', color=color, label=alpha_i)
        axs1[1].plot(J_i_linspace, C_Fy_i_fit, '-', color=color)

        for ax in axs1:
            ax.set_xlabel("$ J $")

        # --- Fig 2: Experiment C_F with quadratic fit lines vs J ---
        ax2.plot(J_i_flat, C_F_i_flat, 'x', color=color, label=alpha_i)
        ax2.plot(J_i_linspace, C_F_i_fit, '-', color=color)

        # --- Model evaluation ---
        phi = 0.6
        sigma = 1.3
        Mb_model = 0.05
        gam = 1.4

        C_Tx_model = phi * (phi / sigma * sin_alpha_i - J_i_linspace)
        C_Ty_model = np.ones_like(J_i_linspace) * phi**2 / sigma * cos_alpha_i
        C_T_model = C_Ty_model * cos_alpha_i + C_Tx_model * sin_alpha_i

        # --- Fig 3: Model C_Tx and C_Ty and experiment C_Fx and C_Fy vs J ---
        axs3[0].plot(J_i_linspace, C_Tx_model, '-', color=color, label=alpha_i)
        axs3[0].plot(J_i_linspace, C_Fx_i_fit, '--', color=color)

        axs3[1].plot(J_i_linspace, C_Ty_model, '-', color=color, label=alpha_i)
        axs3[1].plot(J_i_linspace, C_Fy_i_fit, '--', color=color)

        # --- Fig 4: Model C_T and experiment C_F vs J ---
        ax4.plot(J_i_linspace, C_T_model, '-', color=color, label=alpha_i)
        ax4.plot(J_i_linspace, C_F_i_fit, '--', color=color)

        # --- Fig 5: Difference between model C_Ti and experiment C_Fi ---
        C_Fx_i_trunc_fit = C_Fx_i_fit_func.convert().cutdeg(1)(J_i_linspace)
        C_Fy_i_trunc_fit = C_Fy_i_fit_func.convert().cutdeg(1)(J_i_linspace)

        # axs5[0].plot(J_i_linspace, C_Fx_i_fit, '--', color=color, label=alpha_i)
        axs5[0].plot(J_i_linspace, C_Fx_i_trunc_fit, '--', color=color, label=alpha_i)
        axs5[0].plot(J_i_linspace, C_Tx_model, '-', color=color, label=alpha_i)

        # axs5[1].plot(J_i_linspace, C_Fy_i_fit, '--', color=color, label=alpha_i)
        axs5[1].plot(J_i_linspace, C_Fy_i_trunc_fit, '--', color=color, label=alpha_i)
        axs5[1].plot(J_i_linspace, C_Ty_model, '-', color=color, label=alpha_i)

        # axs5[0].plot(J_i_linspace, C_Fx_i_fit - C_Tx_model, '-', color=color, label=alpha_i)
        # axs5[1].plot(J_i_linspace, C_Fy_i_fit - C_Ty_model, '-', color=color, label=alpha_i)


axs1[0].legend()
axs1[1].legend()
ax2.legend()
axs3[0].legend()
axs3[1].legend()
ax4.legend()
plt.show()
