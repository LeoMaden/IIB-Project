import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
from dataclasses import dataclass




@dataclass
class ExpGeom:
    name : str
    rhub : float
    rtip : float 
    rmean : float
    A1 : float

@dataclass
class ExpAtm:
    Pa : float
    Ta : float
    Ha : float

@dataclass
class ExpRaw:
    V_supply : float
    P : np.ndarray
    T_motor : np.ndarray
    I_supply : np.ndarray
    shaft_freq : np.ndarray
    V_long_avg : np.ndarray
    V_lat_avg : np.ndarray

@dataclass
class ExpData:
    geom : ExpGeom
    atm : ExpAtm
    raw : ExpRaw
    aoa : np.ndarray
    tunnel_switch : np.ndarray
    pwm_speed : np.ndarray
    load_calib : np.ndarray

    def __post_init__(self):
        self.geom = ExpGeom(**self.geom)
        self.atm = ExpAtm(**self.atm)
        self.raw = ExpRaw(**self.raw)


def parse(s):
    names = s.dtype.names

    s_dict = {}

    for name in names:
        content = s[name].item()

        if type(content) == np.ndarray and content.shape == ():
            # Data stored in s.name is a struct; recurse
            s_dict[name] = parse(content)
        else:
            # Content is primitive type of ndarray
            s_dict[name] = content

    return s_dict

def load_data():
    e_runs = {}
    N_runs = {}
    runs = ["run1", "run2"]
    path = "Experimental Data/{}_proc.mat"

    for run in runs:
        run_path = path.format(run)
        data = scipy.io.loadmat(run_path, squeeze_me=True)

        e_runs[run] = parse(data["e"])
        N_runs[run] = parse(data["N"])

    return e_runs, N_runs

if __name__ == "__main__":
    e_runs, N_runs = load_data()

    e = ExpData(**e_runs["run1"])
    pass