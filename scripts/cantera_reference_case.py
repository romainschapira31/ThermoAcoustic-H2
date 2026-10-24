"""
Reference case: laminar flame speed and adiabatic temperature of a
stoichiometric CH4/air mixture (298.15 K, 1 atm), compared with the
compiled measurement of FFCM-1 (36 +/- 1.7 cm/s).

The parameters are fixed on purpose: this script reproduces a reference
case, so they must not be edited.
"""

import cantera as ct

# CONSTANTS
T_in = 298.15
P_in = 101325
fuel = "CH4"
oxidizer = "O2:1, N2:3.76"
width = 0.03
phi = 1.0
refine_criteria = {"ratio": 3, "slope": 0.05, "curve": 0.05}


gas = ct.Solution('gri30.yaml')

gas.TP = T_in, P_in

gas.set_equivalence_ratio(phi, fuel, oxidizer)

flame = ct.FreeFlame(gas, width = width)


flame.set_refine_criteria(**refine_criteria)

flame.solve(loglevel=0, auto=True)

flame_speed = flame.velocity[0]
T_adiabatic = flame.T[-1]

print(f"The flame's speed = {flame_speed:.4f} m/s")
print(f"The maximal temperature = {T_adiabatic:.2f} K ")
