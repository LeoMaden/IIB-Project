import pandas as pd
import numpy as np


def optimum_Cpr2_symmetrical(A2_A1):
    x = A2_A1
    xp = Cpr2_data_symmetrical["A2/A1"]

    fp1 = Cpr2_data_symmetrical["Cpr"]
    fp2 = Cpr2_data_symmetrical["L/h1"]

    Cpr2 = np.interp(x, xp, fp1)
    L_h1 = np.interp(x, xp, fp2)

    return Cpr2, L_h1

def calc_phi(A2_A1, L_h1):
    tan_phi = (A2_A1 - 1) / (2 * L_h1)
    phi = np.rad2deg(np.arctan(tan_phi))
    return phi


# Load data
#Cpr1_data_symmetrical = pd.read_csv("ESDU Data/Cpr*_fig3.csv", header=0)
Cpr2_data_symmetrical = pd.read_csv("ESDU Data/Cpr**_fig3.csv", header=0)