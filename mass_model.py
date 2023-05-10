from dataclasses import dataclass

@dataclass
class MassParameters:
    cowl_density: float
    cowl_solidity: float
    hub_density: float
    hub_solidity: float
    rotor_density: float
    stator_density: float