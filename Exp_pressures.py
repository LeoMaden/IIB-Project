from modules import experimental
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 20})


data = experimental.load_data(["Experimental Data/run3_processed.mat"])

for run in data:
    P_static = run.P_static
    P_total = run.P_total

    R = 287
    rho = run.atm.Pa / (R * run.atm.Ta)

    V = np.sqrt(2 * (run.P_total - run.P_static) / rho)
    U = 2 * np.pi * run.f_shaft * run.geom.rmean
    J = V / U

    gam = 1.4
    M = V / np.sqrt(gam * R * run.atm.Ta)
    Mb = U / np.sqrt(gam * R * run.atm.Ta)
    P1 = run.atm.Pa + P_static
    P1_nd = P1 / (0.5 * rho * V**2)

    V_flat = V.flatten()
    P1_nd_flat = P1_nd.flatten()
    J_flat = J.flatten()

    fig, ax = plt.subplots()

    ax.plot(J_flat, P1_nd_flat, "x")
    
    i_valid = np.isfinite(J_flat) & np.isfinite(P1_nd_flat)

    min_J = min(J_flat[i_valid])
    max_J = max(J_flat[i_valid])
    max_J = 1
    J_model = np.linspace(min_J, max_J, 500)
    #M_model = V_model / np.sqrt(gam * R * run.atm.Ta)
    Mb_model = 0.05
    P1_nd_model = 2 / (gam * J_model**2 * Mb_model**2)

    ax.plot(J_model, P1_nd_model, "-")
    ax.set_xlim(0, 0.1)
    ax.set_ylim(0, 1e6)

    # M = np.linspace(0.01, 0.05, 100)
    # P1_nd = 2 / (gam * M**2)

    # ax.plot(M, P1_nd)


plt.show()