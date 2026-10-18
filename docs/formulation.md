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

The acoustic field inside the duct network is described by 1D linear
acoustics with no mean flow. At any point, the state is the pair (p, u):
acoustic pressure and acoustic velocity. Each duct element maps the state
at its inlet to the state at its outlet via a 2x2 transfer matrix:

    [p]      [ cos(kL)        -i*Z*sin(kL) ] [p]
    [u]_out =[ -i*sin(kL)/Z    cos(kL)     ] [u]_in

where k = omega/c is the wavenumber, L the duct length, and Z = rho*c the
characteristic acoustic impedance of the gas.

Boundary conditions are expressed in the general form a*p + b*u = 0:
- Closed end (rigid wall): u = 0, i.e. (a, b) = (0, 1)
- Open end (idealized, no end correction, no radiation impedance): p = 0,
  i.e. (a, b) = (1, 0)

Eigenfrequencies of the network are found by propagating a state that
satisfies the inlet boundary condition through all elements, and searching
for frequencies at which the outlet boundary condition is also satisfied
(characteristic equation). A frequency sweep is used to locate the minima
of the characteristic value's magnitude.

### Validation: closed-open duct

The solver was validated on a simple closed-open duct (L = 1 m,
c = 343 m/s), for which the analytical eigenfrequencies are known:

f_n = (2n-1) * c / (4L), n = 1, 2, 3, ...

| Mode | Numerical (Hz) | Analytical (Hz) | Relative error (%) |
|---|---|---|---|
| 1 | 86.0 | 85.75 | 0.292 |
| 2 | 257.0 | 257.25 | 0.097 |
| 3 | 429.0 | 428.75 | 0.058 |
| 4 | 600.0 | 600.25 | 0.042 |
| 5 | 772.0 | 771.75 | 0.032 |
| 6 | 943.0 | 943.25 | 0.027 |

<img src="../results/figures/acoustics_validation.png" alt="Acoustic resonance search for a closed-open duct" width="600">

All six modes are recovered within 0.3% of the analytical values, confirming
the correctness of the transfer-matrix implementation. The residual error
is purely numerical (finite resolution of the frequency sweep, 1 Hz step
over 1-1000 Hz) rather than a modeling limitation.

This validates the acoustic-only building block of the network, isolated
from any flame response. The full thermoacoustic model (coupled acoustics + flame response) will be validated once the n-tau flame model is implemented (Section 3).

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
