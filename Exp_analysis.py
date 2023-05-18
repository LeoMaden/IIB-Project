"""
Analyse wind tunnel data and compute nondimensionals
for comparison to CFD and use in the model
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from modules import experimental
from dataclasses import dataclass

@dataclass
class ExpAnalysis:
    rho: float
    U: np.ndarray
    A3: float

    J: np.ndarray
    C_Fx: np.ndarray
    C_Fy: np.ndarray
    C_dotWx: np.ndarray

def analyse_run(run_data: experimental.ExpData):
    R = 287     # Specific gas constant for air
    rho = run_data.atm.Pa / (R * run_data.atm.Ta)
    V = np.sqrt(2 * (run_data.P_total - run_data.P_static) / rho)
    U = run_data.f_shaft * run_data.geom.rmean

    data = {}

    data["U"] = U
    data["rho"] = rho

    # Area at rotor inlet
    A3 = run_data.geom.A1
    data["A3"] = A3

    # Airspeed ratio
    data["J"] = V / U

    # Force components: x is aligned to freestream and
    # y is perpendicular 
    Fx = run_data.F_long
    Fy = run_data.F_lat 

    # Force coefficients defined same way as thrust coefficients
    data["C_Fx"] = Fx / (rho * A3 * U**2)
    data["C_Fy"] = Fy / (rho * A3 * U**2)

    # Electrical power
    dotWx = run_data.V_supply * run_data.I_supply

    # Power coefficent 
    data["C_dotWx"] = dotWx / (rho * A3 * U**3)

    return ExpAnalysis(**data)


# Load runs and analyse data
runs = ["run2"]
file_paths = ["Experimental Data/{}_processed.mat".format(run) for run in runs]

e_runs = experimental.load_data(file_paths)
runs_calcs = [analyse_run(e_run) for e_run in e_runs]

# Plot force coefficients when tunnel is off, in this case
# they corresond to the thrust coefficients as there is no drag
plt.figure(figsize=(10, 8))
markers = ["x", "o", "s", "+", "D", "^", "*"]
lines = ["-", "--"]
colours = mcolors.TABLEAU_COLORS.keys()

for i_run, run_calcs in enumerate(runs_calcs):
    e_run = e_runs[i_run]

    # Find index where tunnel is off
    i_off = np.where(e_run.tunnel_switch == 0)[0]

    alpha = e_run.aoa
    rpm = 60 * e_run.f_shaft

    C_Tx = np.squeeze(run_calcs.C_Fx[:, i_off, :])
    C_Ty = np.squeeze(run_calcs.C_Fy[:, i_off, :])

    # alpha, C_Tx, C_Ty, pwm_speed
    for i in range(len(e_run.pwm_speed)):
        plt.plot(alpha, C_Tx[:, i], marker=markers[i], color="b")
        plt.plot(alpha, C_Ty[:, i], marker=markers[i], color="r")

plt.show()


# if __name__ == "__main__":

#     runs = ["run1", "run2"]
#     file_paths = ["Experimental Data/{}_processed.mat".format(run) for run in runs]

#     e_runs, N_runs = load_data(file_paths)

#     print(N_runs)
# pass