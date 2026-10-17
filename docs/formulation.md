# Physical Formulation

## 1. Chemical kinetics and laminar flame calibration

The laminar flame speed and adiabatic flame temperature of the fresh CH4/air
and H2/air mixtures are computed with Cantera, using a freely-propagating
premixed flame (`FreeFlame`) solved with the GRI-Mech 3.0 mechanism for
methane and a dedicated H2/air mechanism for hydrogen.

All calculations in this project (validation case and main thermoacoustic
study) use T = 298.15 K and p = 1 atm (101325 Pa) as the reference fresh-gas
conditions, matching standard laboratory conditions (25degC) used throughout
the combustion literature.

These quantities are not the end goal in themselves: they are used later to
calibrate the flame response model (n-tau or flame transfer function) that
drives the thermoacoustic network.

### Reference case: CH4/air, phi = 1.0, T = 298.15 K, P = 1 atm

| Quantity | Value |
|---|---|
| Laminar flame speed | 0.38 m/s |
| Adiabatic flame temperature | 2229.5 K |

### Numerical setup

- Domain width: 0.03 m
- Refine criteria: `ratio=3`, `slope=0.05`, `curve=0.05`
- Solved with `auto=True` (staged energy-off / energy-on / refine)

### Grid convergence study

| Grid points | Refine criteria | Flame speed | Adiabatic T |
|---|---|---|---|
| 34 | Cantera default | 0.42 m/s | 2225.8 K |
| 139 | ratio=3, slope=curve=0.1 | 0.38 m/s | 2229.3 K |
| 253 | ratio=3, slope=curve=0.05 | 0.38 m/s | 2229.5 K |

The flame speed is unchanged between the 139-point and 253-point grids,
confirming grid independence. The default Cantera refine criteria (34 points)
under-resolve the flame front and overpredict the flame speed by ~10%.

The adiabatic temperature is far less sensitive to grid resolution than the
flame speed: it is a global energy-balance quantity, nearly independent of
the local spatial discretization, while the flame speed depends directly on
the gradients and reaction rates inside the thin flame front, which require
a finer mesh to capture accurately.

### Validation against reference

The computed flame speed (0.38 m/s) is compared to the experimental value of
Vagelopoulos & Egolfopoulos (1998), 0.36 m/s +/- 0.017 m/s, measured by the
counterflow twin-flame technique at the same conditions (CH4/air, phi=1.0,
298 K, 1 atm):

> Vagelopoulos, C. M., & Egolfopoulos, F. N. (1998). Direct experimental
> determination of laminar flame speeds. Proceedings of the Combustion
> Institute, 27, 513-519.

The ~5.6% overprediction is consistent with GRI-Mech 3.0's known tendency to
slightly overpredict flame speed near stoichiometry, and is within the
typical scatter (+/-2-3 cm/s) observed between different experimental
techniques reported in the literature (Bunsen burner, heat-flux method,
counterflow, spherical/cylindrical constant-volume vessels).

The computed adiabatic flame temperature (2229.5 K) was cross-checked
against an independent Cantera/GRI-Mech 3.0 implementation (agreement within
0.2%), confirming the correctness of this implementation:

> Adiabatic flame temperatures for oxy-methane, oxy-hydrogen, air-methane,
> and air-hydrogen stoichiometric combustion using the NASA CEARUN tool,
> GRI-Mech 3.0 reaction mechanism, and Cantera Python package. arXiv:2503.11826.

As a physical sanity check, this result is also ~100 K lower than a
simplified complete-combustion estimate without dissociation (2325.6 K,
based on Turns' textbook thermodynamic tables). This is expected: at flame
temperature, dissociation of combustion products (CO, OH, H, O) absorbs
energy and lowers the true equilibrium temperature below the simplified
no-dissociation estimate.

## 2. Acoustic model

(1D linear acoustics, transfer matrix approach — to complete once
`acoustics.py` is implemented)

## 3. Flame response model

(n-tau model, Crocco's law — to complete)

## 4. Stability criterion

(Rayleigh criterion — to complete)


## References

1. Vagelopoulos, C. M., & Egolfopoulos, F. N. (1998). Direct experimental
   determination of laminar flame speeds. *Proceedings of the Combustion
   Institute*, 27, 513-519.
   — Experimental reference for the CH4/air laminar flame speed
   (counterflow twin-flame technique).

2. Adiabatic flame temperatures for oxy-methane, oxy-hydrogen, air-methane,
   and air-hydrogen stoichiometric combustion using the NASA CEARUN tool,
   GRI-Mech 3.0 reaction mechanism, and Cantera Python package.
   arXiv:2503.11826.
   — Independent Cantera/GRI-Mech 3.0 implementation used to cross-check
   the correctness of this project's flame calibration script.

3. Turns, S. R. *An Introduction to Combustion: Concepts and Applications.*
   — Simplified complete-combustion (no dissociation) adiabatic flame
   temperature estimate, used as a physical sanity check.
