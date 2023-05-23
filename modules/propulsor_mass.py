from dataclasses import dataclass
from . import geometry


@dataclass
class NonDimensionalMassParams:
    cowl_density: float
    hub_density: float
    rotor_density: float
    stator_density: float
    cowl_solidity: float
    hub_solidity: float

@dataclass
class NonDimensionalMassModel:

    mass_params: NonDimensionalMassParams
    geom: geometry.NonDimensionalGeometry

    def calc_M_rotor(self):
        return self.mass_params.rotor_density * self.geom.calc_Vh3_rotor() * self.geom.A3D_hr3**(-1)
    
    def calc_M_stator(self):
        return self.mass_params.stator_density * self.geom.calc_Vh3_stator() * self.geom.A3D_hr3**(-1)
    
    def calc_M_hub(self, xr_hub):
        Vh3_hub = self.geom.calc_Vh3_hub(xr_hub)
        return self.mass_params.hub_solidity * self.mass_params.hub_density * Vh3_hub * self.geom.A3D_hr3**(-1)
    
    def calc_M_cowl(self, xr_cas, xr_cowl):
        Vh3_cowl = self.geom.calc_Vh3_cowl(xr_cas, xr_cowl)
        return self.mass_params.cowl_solidity * self.mass_params.cowl_density * Vh3_cowl * self.geom.A3D_hr3**(-1)
    

    def calc_W_bat(self):
        pass