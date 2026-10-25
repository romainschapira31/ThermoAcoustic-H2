"""
Reference check: sensitivity of the laminar flame speed to numerical and
transport settings, for CH4/air and H2/air at phi = 1.0, 298.15 K, 1 atm.

It explains why the computed hydrogen flame speed (Cantera's default
mixture-averaged transport) is above the compiled measurement, and shows
that the numerical settings, the size of the mechanism and the 0.15 K
difference with the compilations play no role. The parameters are fixed on
purpose: this script reproduces a reference comparison.
"""

import cantera as ct

CASES = {
    "CH4/air": dict(mechanism="gri30.yaml", fuel="CH4", measured=(0.360, 0.017)),
    "H2/air": dict(mechanism="h2o2.yaml", fuel="H2", measured=(2.1117, 0.1099)),
}

VARIANTS = [
    ("Project calculation (mixture-averaged transport)", {}),
    ("Stricter refinement criteria (slope = curve = 0.03)", dict(slope=0.03, curve=0.03)),
    ("Wider domain (0.10 m)", dict(width=0.10)),
    ("Multicomponent transport", dict(transport="multicomponent")),
    ("Multicomponent transport + thermal diffusion (Soret)", dict(transport="multicomponent", soret=True)),
    ("Full GRI-Mech 3.0 (gri30.yaml)", dict(mechanism="gri30.yaml")),
    ("Fresh gas at 298.0 K (as in the compilations)", dict(temperature=298.0)),
]


def flame_speed(mechanism, fuel, temperature=298.15, width=0.03, ratio=3, slope=0.05,
                curve=0.05, transport=None, soret=False):
    """Laminar flame speed [m/s] of a stoichiometric fuel/air mixture."""
    gas = ct.Solution(mechanism)
    if transport:
        gas.transport_model = transport
    gas.TP = temperature, ct.one_atm
    gas.set_equivalence_ratio(1.0, fuel, "O2:1, N2:3.76")
    flame = ct.FreeFlame(gas, width=width)
    flame.set_refine_criteria(ratio=ratio, slope=slope, curve=curve)
    flame.soret_enabled = soret
    flame.solve(loglevel=0, auto=True)
    return flame.velocity[0]


results = {
    name: [flame_speed(**{"mechanism": case["mechanism"], "fuel": case["fuel"], **settings})
           for _, settings in VARIANTS]
    for name, case in CASES.items()
}

print(f"{'Setting':56s} {'CH4/air':>9s} {'H2/air':>9s}   [m/s]")
for i, (label, _) in enumerate(VARIANTS):
    print(f"{label:56s} {results['CH4/air'][i]:9.3f} {results['H2/air'][i]:9.3f}")
print(f"{'Compiled measurement':56s} "
      f"{CASES['CH4/air']['measured'][0]:9.3f} {CASES['H2/air']['measured'][0]:9.3f}")

print("\nDifference with the compiled measurement [%], in units of its uncertainty:")
for name, case in CASES.items():
    value, uncertainty = case["measured"]
    for i, (label, _) in enumerate(VARIANTS):
        s = results[name][i]
        print(f"  {name:8s} {label:56s} {100 * (s - value) / value:+5.1f} %  "
              f"({(s - value) / uncertainty:+.1f} x)")