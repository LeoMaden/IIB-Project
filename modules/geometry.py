import bezier
from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np

@dataclass
class NonDimensionalGeometry:
    Lh_diffuser:    float   # Diffuser length/rotor height
    Lh_rotor:       float   # Rotor length/rotor height
    Lh_intake:      float   # Intake length/rotor height
    rr_hub_tip:     float   # Hub to tip ratio 
    th_cowl:        float   # Cowl thickness to rotor height ratio
    tc_rotor:       float   # Rotor thickness to chord ratio
    tc_stator:      float   # Stator thickness to chord ratio
    phi_i:          float   # Diffuser inner wall angle
    phi_o:          float   # Diffuser outer wall angle
    cL_rotor:       float
    cL_stator:      float
    delta_oh:       float   # Ratio of delta_o to rotor height
    N_rotor:        float   # Number of rotor blades
    N_stator:       float   # Number of stator blades

    @property
    def Rh_hub(self):
        return (1 / self.rr_hub_tip - 1)**(-1)
    
    @property
    def Rh_cas(self):
        return (1 - self.rr_hub_tip)**(-1)
    
    @property
    def Rh_mean(self):
        return (self.Rh_hub + self.Rh_cas) / 2
    
    @property
    def delta_ih(self):
        return self.Rh_hub - self.Lh_diffuser * np.tan(np.deg2rad(self.phi_i))
    
    @property
    def Lh_total(self):
        return self.Lh_intake + self.Lh_rotor + self.Lh_diffuser
    
    @property
    def A3D_hr3(self):
        return 2 * np.pi * (1 + self.rr_hub_tip) / (1 - self.rr_hub_tip)**2
    
    def calc_Vh3_rotor(self):
        return self.tc_rotor * self.cL_rotor**2 * self.Lh_rotor**2 * self.N_rotor
    
    def calc_Vh3_stator(self):
        tan_phi_i = np.tan(np.deg2rad(self.phi_i))
        tan_phi_o = np.tan(np.deg2rad(self.phi_o))
        x = 1 + self.Lh_diffuser * (tan_phi_i + tan_phi_o)
        return self.tc_stator * self.cL_stator**2 * self.Lh_diffuser**2 * x * self.N_stator

    def calc_Vh3_hub(self, xr_hub):
        x = xr_hub[0, :]
        r = xr_hub[1, :]
        return np.pi * np.trapz(r**2, x)
    
    def calc_Vh3_cowl(self, xr_cas, xr_cowl):
        x = xr_cas[0, :]
        r_cas = xr_cas[1, :]
        r_cowl = xr_cowl[1, :]
        return np.pi * np.trapz(r_cowl**2 - r_cas**2, x)
    
    def calc_hub_line(self, Ns: int=20):
        Rhub = self.Rh_hub
        Lintake = self.Lh_intake
        Lrotor = self.Lh_rotor
        Ldiff = self.Lh_diffuser
        delta_i = self.delta_ih

        # Intake bezier
        nodes = np.array([
            [0,     0,          0.6*Lintake,    1.0*Lintake ],
            [0,     0.8*Rhub,   1.0*Rhub,       1.0*Rhub    ]
        ])
        intake_curve = bezier.Curve.from_nodes(nodes)

        # Rotor line
        nodes = np.array([
            [0,         1*Lrotor    ],
            [1*Rhub,    1*Rhub      ]
        ])
        nodes[0, :] += Lintake
        rotor_curve = bezier.Curve.from_nodes(nodes)

        # Diffuser line
        nodes = np.array([
            [0,         0.4*Ldiff,  0.6*Ldiff,      1*Ldiff     ],
            [1*Rhub,    1.0*Rhub,   1.0*delta_i,    1*delta_i   ]
        ])
        nodes[0, :] += Lintake + Lrotor
        diffuser_curve = bezier.Curve.from_nodes(nodes)

        splines = [intake_curve, rotor_curve, diffuser_curve]
        s_vals = np.linspace(0, 1, Ns)
        xr_hub = np.empty((2, 0))

        for spline in splines:
            xr = spline.evaluate_multi(s_vals)
            xr_hub = np.hstack([xr_hub, xr])

        return xr_hub
    
    def calc_cas_line(self, Ns: int=20):
        Rcas = self.Rh_cas
        tcowl = self.th_cowl
        Lintake = self.Lh_intake
        Lrotor = self.Lh_rotor
        Ldiff = self.Lh_diffuser

        # Intake bezier
        nodes = np.array([
            [0,         0,          0.4*Lintake,    1*Lintake   ],
            [0.3*tcowl, 0.1*tcowl,  0,              0           ]
        ])
        intake_curve = bezier.Curve.from_nodes(nodes)

        # Rotor line
        nodes = np.array([
            [0,     1*Lrotor    ],
            [0,     0           ]
        ])
        nodes[0, :] += Lintake
        rotor_curve = bezier.Curve.from_nodes(nodes)

        # Diffuser line
        y = np.tan(np.deg2rad(self.phi_o)) * Ldiff
        nodes = np.array([
            [0,     0.4*Ldiff,  0.6*Ldiff,  1*Ldiff ],
            [0,     0,          1*y,        1*y     ]
        ])
        nodes[0, :] += Lintake + Lrotor
        diffuser_curve = bezier.Curve.from_nodes(nodes)

        splines = [intake_curve, rotor_curve, diffuser_curve]
        s_vals = np.linspace(0, 1, Ns)
        xr_cas = np.empty((2, 0))

        for spline in splines:
            xr = spline.evaluate_multi(s_vals)
            xr_cas = np.hstack([xr_cas, xr])

        # Translate up
        xr_cas[1, :] += Rcas

        return xr_cas
    
    def calc_cowl_line(self, Ns: int=30):
        Rcas = self.Rh_cas
        tcowl = self.th_cowl
        Lintake = self.Lh_intake
        Ldiff = self.Lh_diffuser
        Ltot = self.Lh_total
        delta_o = self.delta_oh

        # Front line
        nodes = np.array([
            [0,         0,          0.2*Ltot,   0.5*Ltot    ],
            [0.3*tcowl, 0.8*tcowl,  1*tcowl,    1*tcowl     ]
        ])
        front_curve = bezier.Curve.from_nodes(nodes)

        # Rear line
        y = np.tan(np.deg2rad(self.phi_o)) * Ldiff
        nodes = np.array([
            [0,       0.2*Ltot,     0.4*Ltot,   0.5*Ltot        ],
            [1*tcowl, 1*tcowl,      0.5*tcowl,  y + delta_o     ]
        ])
        nodes[0, :] += 0.5*Ltot
        rear_curve = bezier.Curve.from_nodes(nodes)

        splines = [front_curve, rear_curve]
        s_vals = np.linspace(0, 1, Ns)
        xr_cowl = np.empty((2, 0))

        for spline in splines:
            xr = spline.evaluate_multi(s_vals)
            xr_cowl = np.hstack([xr_cowl, xr])

        # Translate up
        xr_cowl[1, :] += Rcas

        return xr_cowl

    # Example nondimensional geometry
    @staticmethod
    def example():
        return NonDimensionalGeometry(
            Lh_diffuser=0.3,    
            Lh_rotor=0.2,     
            Lh_intake=0.1,   
            rr_hub_tip=0.15, 
            th_cowl=0.3,
            tc_rotor=0.1,
            tc_stator=0.1,
            phi_i=20,
            phi_o=20,
            cL_rotor=0.9,
            cL_stator=0.6,
            delta_oh=0.02,
            N_rotor=6,
            N_stator=7
        )

# @dataclass
# class Geometry:
#     L_intake: float
#     L_rotor: float
#     L_diffuser: float
#     R_hub: float
#     R_cas: float
#     t_cowl: float
#     phi_i: float
#     phi_o: float
#     c_rotor: float
#     c_stator: float
#     s_rotor: float
#     s_stator: float
#     t_rotor: float
#     t_stator: float
#     delta_o: float
#     nd_geom: NonDimensionalGeometry
#     D: float

#     # Calculated geometry values
#     @property
#     def R_mean(self):
#         return (self.R_hub + self.R_cas) / 2
    
#     @property
#     def h_rotor(self):
#         return self.R_cas - self.R_hub
    
#     @property
#     def h_exit(self):
#         a = self.L_diffuser * np.tan(np.deg2rad(self.phi_o))
#         b = self.L_diffuser * np.tan(np.deg2rad(self.phi_i))
#         return self.h_rotor + a + b
    
#     @property
#     def h_stator(self):
#         return (self.h_rotor + self.h_exit) / 2
    
#     @property
#     def delta_i(self):
#         return self.R_hub - self.L_diffuser * np.tan(np.deg2rad(self.phi_i))
    
#     @property
#     def L_total(self):
#         return self.L_intake + self.L_rotor + self.L_diffuser
    
#     @property
#     def N_rotor(self):
#         circum = 2 * np.pi * self.R_mean
#         N = circum / self.s_rotor
#         return round(N)
    
#     @property
#     def N_stator(self):
#         circum = 2 * np.pi * self.R_mean
#         N = circum / self.s_stator
#         return round(N)
    
#     def __init__(self, nd_geom: NonDimensionalGeometry, D:float):
#         self.D = D
#         self.nd_geom = nd_geom
#         self._calc_geometry(nd_geom, D)

    
#     def calc_hub_line(self, Ns : int = 20):
#         xr_hub = self.nd_geom.calc_hub_line(Ns)
#         xr_hub *= self.h_rotor
#         return xr_hub
    
#     def calc_cas_line(self, Ns : int = 20):
#         xr_cas = self.nd_geom.calc_cas_line(Ns)
#         xr_cas *= self.h_rotor
#         return xr_cas
    
#     def calc_cowl_line(self, Ns : int = 20):
#         xr_cowl = self.nd_geom.calc_cowl_line(Ns)
#         xr_cowl *= self.h_rotor
#         return xr_cowl

#     # Calculate geometry class from nondimensional geometry
#     def _calc_geometry(self, nd_geom: NonDimensionalGeometry, D: float):
#         """Calculate the dimensional geometry parameters from the non-dimensional
#         parameters and the diameter and set fields 

#         Args:
#             nd_geom (NonDimensionalGeometry): The class containing the values of the non-dimensional geometry

#             D (float): The rotor tip diameter
#         """

#         # Radii
#         self.R_cas = D / 2
#         self.R_hub = nd_geom.rr_hub_tip * self.R_cas

#         # Passage height
#         h_rotor = self.R_cas - self.R_hub

#         # Section lengths
#         self.L_intake = nd_geom.Lh_intake * h_rotor
#         self.L_rotor = nd_geom.Lh_rotor * h_rotor
#         self.L_diffuser = nd_geom.Lh_diffuser * h_rotor

#         # Cowl thickness
#         self.t_cowl = nd_geom.th_cowl * h_rotor

#         # Diffuser angles
#         self.phi_i = nd_geom.phi_i
#         self.phi_o = nd_geom.phi_o

#         # Chords
#         self.c_rotor = nd_geom.cL_rotor * self.L_rotor
#         self.c_stator = nd_geom.cL_stator * self.L_diffuser

#         # Pitch at mean raduis
#         self.s_rotor = nd_geom.sc_rotor * self.c_rotor
#         self.s_stator = nd_geom.sc_stator * self.c_stator

#         # Blade thickness
#         self.t_rotor = nd_geom.tc_rotor * self.c_rotor
#         self.t_stator = nd_geom.tc_stator * self.c_stator

#         # Delta o
#         self.delta_o = nd_geom.delta_oh * h_rotor
    
#     @staticmethod
#     def example():
#         D = 0.5
#         nd_geom = NonDimensionalGeometry.example()
#         geom = Geometry(nd_geom, D)
#         return geom
