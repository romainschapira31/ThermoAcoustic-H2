"""
Grid convergence of the laminar flame speed and adiabatic temperature of a
stoichiometric CH4/air flame (GRI-Mech 3.0, 298.15 K, 1 atm).

The three refinement settings are the Cantera default and the two stricter
criteria used by the reference case (cantera_reference_case.py). The
parameters are fixed on purpose: this script reproduces a reference study.
"""

import cantera as ct

SETTINGS = [
    ("Cantera default", None),
    ("ratio=3, slope=curve=0.1", dict(ratio=3, slope=0.1, curve=0.1)),
    ("ratio=3, slope=curve=0.05", dict(ratio=3, slope=0.05, curve=0.05)),
]

print(f"{'Refinement criteria':28s} {'points':>7s} {'S_L [m/s]':>10s} {'T_ad [K]':>10s}")
for label, criteria in SETTINGS:
    gas = ct.Solution("gri30.yaml")
    gas.TP = 298.15, ct.one_atm
    gas.set_equivalence_ratio(1.0, "CH4", "O2:1, N2:3.76")
    flame = ct.FreeFlame(gas, width=0.03)
    if criteria:
        flame.set_refine_criteria(**criteria)
    flame.solve(loglevel=0, auto=True)
    print(f"{label:28s} {len(flame.grid):7d} {flame.velocity[0]:10.4f} {flame.T[-1]:10.1f}")