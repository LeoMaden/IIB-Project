import pandas as pd
import numpy as np


def optimum_Cpr2_symmetrical(A2_A1):
    # Add point at A2/A1 = 1, Cpr = 0, L/h1 = 0
    new_data = pd.DataFrame.from_dict({"A2/A1": [1], "Cpr": [0], "L/h1": [0]})
    data = pd.concat([Cpr2_data_symmetrical, new_data]).sort_values("A2/A1")

    x = A2_A1
    xp = data["A2/A1"]

    fp1 = data["Cpr"]
    fp2 = data["L/h1"]

    Cpr2 = np.interp(x, xp, fp1, left=0)
    L_h1 = np.interp(x, xp, fp2, left=0)

    return Cpr2, L_h1

def calc_phi(A2_A1, L_h1):
    if L_h1 == 0:
        return 0

    tan_phi = (A2_A1 - 1) / (2 * L_h1)
    phi = np.rad2deg(np.arctan(tan_phi))
    return phi


# Load data
#Cpr1_data_symmetrical = pd.read_csv("ESDU Data/Cpr*_fig3.csv", header=0)
Cpr2_data_symmetrical = pd.read_csv("ESDU Data/Cpr**_fig3.csv", header=0)