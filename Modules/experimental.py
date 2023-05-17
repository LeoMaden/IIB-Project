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

# @dataclass
# class ExpRaw:
#     V_supply : float
#     P : np.ndarray
#     T_motor : np.ndarray
#     I_supply : np.ndarray
#     shaft_freq : np.ndarray
#     V_long_avg : np.ndarray
#     V_lat_avg : np.ndarray

# @dataclass
# class ExpData:
#     geom : ExpGeom
#     atm : ExpAtm
#     raw : ExpRaw
#     aoa : np.ndarray
#     tunnel_switch : np.ndarray
#     pwm_speed : np.ndarray
#     load_calib : np.ndarray

#     def __post_init__(self):
#         self.geom = ExpGeom(**self.geom)
#         self.atm = ExpAtm(**self.atm)
#         self.raw = ExpRaw(**self.raw)

@dataclass
class ExpData:
    geom : ExpGeom
    atm : ExpAtm
    
    V_supply : float
    I_supply : np.ndarray
    P_total : np.ndarray
    P_static : np.ndarray
    T_motor : np.ndarray
    f_shaft : np.ndarray
    F_long : np.ndarray
    F_lat : np.ndarray

    aoa : np.ndarray
    tunnel_switch : np.ndarray
    pwm_speed : np.ndarray

    def __post_init__(self):
        self.geom = ExpGeom(**self.geom)
        self.atm = ExpAtm(**self.atm)

    @staticmethod
    def from_dict(data):
        return ExpData(**data)



def parse_mat(data):
    names = data.dtype.names

    s_dict = {}

    for name in names:
        content = data[name].item()

        if type(content) == np.ndarray and content.shape == ():
            # Data stored in s.name is a struct; recurse
            s_dict[name] = parse_mat(content)
        else:
            # Content is primitive type or non empty ndarray
            s_dict[name] = content

    return s_dict

def load_data(file_paths):
    e_runs = []

    for i, path in enumerate(file_paths):
        data = scipy.io.loadmat(path, squeeze_me=True)        
        
        e_dict = parse_mat(data["e"])
        exp_data = ExpData.from_dict(e_dict)
        e_runs.append(exp_data)

    return e_runs