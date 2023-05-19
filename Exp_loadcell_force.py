from modules import experimental
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})


data = experimental.load_data(["Experimental Data/run3_processed.mat"])

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 8), sharey="all", sharex="all")
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
    
    i = 4
    x = J[i, :, :].squeeze().T
    y1 = C_Fx[i, :, :].squeeze().T
    y2 = C_Fy[i, :, :].squeeze().T


    fig.suptitle(fr"$ \alpha = $ {alpha[i]}")

    axs[0].plot(x, y1, "x")
    axs[0].set_title("$ C_{F_x} $")

    axs[1].plot(x, y2, "x")
    axs[1].set_title("$ C_{F_y} $")

    for ax in axs:
        ax.legend(run.tunnel_switch)
        ax.set_xlabel("$ J $")

plt.show()
