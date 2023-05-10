import bezier
from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np

@dataclass
class NonDimensionalGeometry:
    Lh_diffuser:    float   # Diffuser length/rotor height
    Lh_rotor:       float   # Rotor length/rotor height
    Lh_intake:      float   # Intake length/rotor height
    sc_rotor:       float   # Rotor pitch to chord ratio
    sc_stator:      float   # Stator pitch to chord ratio
    rr_hub_tip:     float   # Hub to tip ratio 
    th_cowl:        float   # Cowl thickness to rotor height ratio
    tc_rotor:       float   # Rotor thickness to chord ratio
    tc_stator:      float   # Stator thickness to chord ratio
    phi_i:          float   # Diffuser inner wall angle
    phi_o:          float   # Diffuser outer wall angle
    cL_rotor:       float
    cL_stator:      float
    delta_oh:       float   # Ratio of delta_o to rotor height

    # Example nondimensional geometry
    @staticmethod
    def example():
        return NonDimensionalGeometry(
            Lh_diffuser=0.3,    
            Lh_rotor=0.2,     
            Lh_intake=0.1,   
            sc_rotor=5,  
            sc_stator=5,     
            rr_hub_tip=0.15, 
            th_cowl=0.3,
            tc_rotor=0.1,
            tc_stator=0.1,
            phi_i=20,
            phi_o=20,
            cL_rotor=0.9,
            cL_stator=0.6,
            delta_oh=0.02
        )

@dataclass
class Geometry:
    L_intake: float
    L_rotor: float
    L_diffuser: float
    R_hub: float
    R_cas: float
    t_cowl: float
    phi_i: float
    phi_o: float
    c_rotor: float
    c_stator: float
    s_rotor: float
    s_stator: float
    t_rotor: float
    t_stator: float
    delta_o: float

    # Calculated geometry values
    @property
    def R_mean(self):
        return (self.R_hub + self.R_cas) / 2
    
    @property
    def h_rotor(self):
        return self.R_cas - self.R_hub
    
    @property
    def delta_i(self):
        return self.R_hub - self.L_diffuser * np.tan(np.deg2rad(self.phi_i))
    
    @property
    def L_total(self):
        return self.L_intake + self.L_rotor + self.L_diffuser
    
    @property
    def N_rotor(self):
        circum = 2 * np.pi * self.R_mean
        N = circum / self.s_rotor
        return round(N)
    
    @property
    def N_stator(self):
        circum = 2 * np.pi * self.R_mean
        N = circum / self.s_stator
        return round(N)

    # Calculate geometry class from nondimensional geometry
    @staticmethod
    def calc_geometry(nd_geom: NonDimensionalGeometry, D: float):
        """Calculate the dimensional geometry parameters from the non-dimensional
        parameters and the diameter

        Args:
            nd_geom (NonDimensionalGeometry): The class containing the values of the non-dimensional geometry

            D (float): The rotor tip diameter
        """

        # Radii
        R_cas = D / 2
        R_hub = nd_geom.rr_hub_tip * R_cas

        # Passage height
        h_rotor = R_cas - R_hub

        # Section lengths
        L_intake = nd_geom.Lh_intake * h_rotor
        L_rotor = nd_geom.Lh_rotor * h_rotor
        L_diffuser = nd_geom.Lh_diffuser * h_rotor

        # Cowl thickness
        t_cowl = nd_geom.th_cowl * h_rotor

        # Diffuser angles
        phi_i = nd_geom.phi_i
        phi_o = nd_geom.phi_o

        # Chords
        c_rotor = nd_geom.cL_rotor * L_rotor
        c_stator = nd_geom.cL_stator * L_diffuser

        # Pitch at mean raduis
        s_rotor = nd_geom.sc_rotor * c_rotor
        s_stator = nd_geom.sc_stator * c_stator

        # Blade thickness
        t_rotor = nd_geom.tc_rotor * c_rotor
        t_stator = nd_geom.tc_stator * c_stator

        # Delta o
        delta_o = nd_geom.delta_oh * h_rotor

        geom = Geometry(
            L_intake=L_intake,
            L_rotor=L_rotor,
            L_diffuser=L_diffuser,
            R_hub=R_hub,
            R_cas=R_cas,
            t_cowl=t_cowl,
            phi_i=phi_i,
            phi_o=phi_o,
            c_rotor=c_rotor,
            c_stator=c_stator,
            s_rotor=s_rotor,
            s_stator=s_stator,
            t_rotor=t_rotor,
            t_stator=t_stator,
            delta_o=delta_o
        )

        return geom
    
    @staticmethod
    def example():
        D = 0.5
        nd_geom = NonDimensionalGeometry.example()
        geom = Geometry.calc_geometry(nd_geom, D)
        return geom


def calc_hub_line(geom: Geometry, Ns: int=20):

    # s values to evaluate beziers at
    s_vals = np.linspace(0, 1, Ns)

    # Intake bezier
    nodes = np.array([
        [0, 0, 0.6, 1],
        [0, 0.8, 1, 1]
    ])
    intake_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([geom.L_intake, geom.R_hub])
    xr_intake = scale @ intake_curve_nd.evaluate_multi(s_vals)

    # Rotor line
    nodes = np.array([
        [0, 1],
        [1, 1]
    ])
    rotor_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([geom.L_rotor, geom.R_hub])
    xr_rotor = scale @ rotor_curve_nd.evaluate_multi(s_vals)
    xr_rotor[0, :] += geom.L_intake

    # Diffuser line
    r = geom.delta_i / geom.R_hub
    nodes = np.array([
        [0, 0.4, 0.6, 1],
        [1, 1, r, r]
    ])
    diffuser_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([geom.L_diffuser, geom.R_hub])
    xr_diffuser = scale @ diffuser_curve_nd.evaluate_multi(s_vals)
    xr_diffuser[0, :] += geom.L_intake + geom.L_rotor

    xr_hub = np.hstack([xr_intake, xr_rotor, xr_diffuser])
    return xr_hub

def calc_cas_line(geom: Geometry, Ns: int=20):

    # s values to evaluate beziers at
    s_vals = np.linspace(0, 1, Ns)

    # Intake bezier
    nodes = np.array([
        [0, 0, 0.4, 1],
        [1, 0.4, 0, 0]
    ])
    intake_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([geom.L_intake, 0.2*geom.t_cowl])
    xr_intake = scale @ intake_curve_nd.evaluate_multi(s_vals)
    xr_intake[1, :] += geom.R_cas

    # Rotor line
    nodes = np.array([
        [0, 1],
        [1, 1]
    ])
    rotor_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([geom.L_rotor, geom.R_cas])
    xr_rotor = scale @ rotor_curve_nd.evaluate_multi(s_vals)
    xr_rotor[0, :] += geom.L_intake

    # Diffuser line
    nodes = np.array([
        [0, 0.4, 0.6, 1],
        [0, 0, 1, 1]
    ])
    diffuser_curve_nd = bezier.Curve.from_nodes(nodes)

    tan_phi_o = np.tan(np.deg2rad(geom.phi_o))
    scale = np.diag([geom.L_diffuser, geom.L_diffuser * tan_phi_o])
    xr_diffuser = scale @ diffuser_curve_nd.evaluate_multi(s_vals)
    xr_diffuser[0, :] += geom.L_intake + geom.L_rotor
    xr_diffuser[1, :] += geom.R_cas

    xr_cas = np.hstack([xr_intake, xr_rotor, xr_diffuser])
    return xr_cas

def calc_cowl_line(geom: Geometry, Ns: int=20):

    # s values to evaluate beziers at
    s_vals = np.linspace(0, 1, Ns)

    # Front bezier
    nodes = np.array([
        [0, 0, 0.4, 1],
        [0, 0.7, 1, 1]
    ])
    front_curve_nd = bezier.Curve.from_nodes(nodes)

    scale = np.diag([0.5*geom.L_total, (1-0.2)*geom.t_cowl])
    xr_front = scale @ front_curve_nd.evaluate_multi(s_vals)
    xr_front[1, :] += geom.R_cas + 0.2*geom.t_cowl

    # Rear bezier
    nodes = np.array([
        [0, 0.4, 0.8, 1],
        [1, 1, 0.4, 0]
    ])
    rear_curve_nd = bezier.Curve.from_nodes(nodes)

    tan_phi_o = np.tan(np.deg2rad(geom.phi_o))
    offset = geom.L_diffuser * tan_phi_o + geom.delta_o
    scale = np.diag([0.5*geom.L_total, geom.t_cowl - offset])

    xr_rear = scale @ rear_curve_nd.evaluate_multi(s_vals)
    xr_rear[0, :] += 0.5*geom.L_total
    xr_rear[1, :] += geom.R_cas + offset

    xr_cowl = np.hstack([xr_front, xr_rear])
    return xr_cowl