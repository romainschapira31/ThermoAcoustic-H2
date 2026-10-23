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
from any flame response. The full thermoacoustic model (coupled acoustics + flame response) is validated in Section 3.

## 3. Flame response model and thermoacoustic coupling

### 3.1 Crocco's n-tau model

The flame's heat release response to acoustic perturbations is modeled
using Crocco's n-tau law: the fluctuating heat release Q' follows the
upstream acoustic velocity u' with a time delay tau and a gain n:

    Q'(t) = n * u'(t - tau)

In the frequency domain, this becomes:

    Q'(omega) = n * exp(-i*omega*tau) * u'(omega)

### 3.2 Compact flame element

The flame is modeled as compact (zero-thickness) within the acoustic
network. Pressure is continuous across the flame, while the acoustic
velocity jumps in proportion to the fluctuating heat release:

    p_out = p_in
    u_out = u_in + coupling_factor * n * exp(-i*omega*tau) * u_in

which gives the 2x2 transfer matrix:

    [p]      [ 1                                          0 ] [p]
    [u]_out =[ 0    1 + coupling_factor * n * exp(-i*omega*tau) ] [u]_in

Eigenfrequencies of the coupled system are now complex: omega = omega_r +
i*omega_i. The growth rate is defined as sigma = -omega_i, so that
sigma > 0 indicates an unstable (growing) mode and sigma < 0 a stable
(decaying) mode, consistent with the Rayleigh criterion.

Roots of the characteristic equation are found with a Newton-Raphson
solver (`scipy.optimize.newton`), starting from an initial guess close to
the acoustic-only eigenfrequency and updated progressively as the flame
parameters are varied.

### 3.3 Consistency check: n = 0 limit

As a sanity check, setting n = 0 (flame heat release does not respond to
acoustic perturbations) should recover the acoustic-only eigenfrequency
found in Section 2, with zero growth rate. This was confirmed numerically:

| Case | omega (rad/s) | Frequency (Hz) | Growth rate |
|---|---|---|---|
| Acoustic-only (Section 2) | 538.78 | 85.75 | 0 (by construction) |
| Coupled, n=0 | 538.78 - 3.8e-21j | 85.75 | ~0 (numerical noise) |

This confirms that the CompactFlame element reduces exactly to a neutral
(non-coupling) element when the flame's heat release response is switched
off, validating the implementation.

### 3.4 Eigenvalue tracking: a methodological caveat

The characteristic equation of the coupled system is transcendental: it
admits infinitely many roots in the complex plane, one per acoustic mode,
further multiplied by the exp(-i*omega*tau) delay term. A Newton solver
returns whichever root lies in the basin of attraction of its initial
guess, with no guarantee that successive solves stay on the same physical
branch.

An initial attempt to build a 2D stability map over the (n, tau) plane
failed for this reason. Sweeping the grid with nearest-neighbour
continuation produced a map that was smooth in appearance and converged at
every point, but was in fact a mosaic of different modes: beyond n ~ 0.9
the solver had silently switched to a spurious branch around 45 Hz,
unrelated to the fundamental mode of the underlying acoustic system
(85.75 Hz).

The discrepancy was detected by cross-checking the map against a direct
1D sweep along n at fixed tau:

| n | Map: f [Hz] | Map: sigma [1/s] | Direct: sigma [1/s] |
|---|---|---|---|
| 0.10 | 87.62 | -5.09 | -5.08 |
| 0.51 | 98.08 | -6.25 | -7.41 |
| 0.93 | (rejected) | — | +35.88 |
| 1.96 | 47.06 | -174.72 | +115.67 |
| 3.00 | 41.64 | -144.60 | +167.54 |

The two methods agree to within rounding at low n, then diverge entirely.
The frequency column identifies the cause unambiguously: the map drifts
from 87 Hz to 45 Hz, i.e. onto a different branch.

The 2D map was therefore discarded. All results below are obtained by
direct 1D sweeps with continuation along a single parameter, and the real
part of omega is monitored throughout as a branch-tracking check.

A rigorous 2D map would require a method that does not rely on local
root-finding — for instance contour-integral eigenvalue methods, or
recasting the problem in matrix eigenvalue form. This is outside the scope
of a low-order study.

### 3.5 Flame calibration from chemistry

The time delay tau is estimated as a convective transit time between the
perturbation location and the flame front:

    tau = d / S_L

where S_L is the laminar flame speed computed with Cantera (Section 1) and
d is a characteristic convective length scale. With d = 10 mm:

| Fuel | phi | S_L [m/s] | tau [ms] |
|---|---|---|---|
| CH4/air | 1.0 | 0.38 | 26.64 |
| H2/air | 0.9 | 2.00 | 4.99 |

The hydrogen flame, being roughly five times faster, yields a delay five
times shorter. This is the only quantity in the flame model that is
derived from chemistry; the interaction index n is not calibrated and is
treated as a free parameter below.

### 3.6 Fuel comparison: growth rate versus interaction index

At each fuel's calibrated tau, the growth rate of the fundamental mode was
computed over n in [0.1, 3.0] by direct root-finding with continuation.

<img src="../results/figures/fuel_comparison_sigma_vs_n.png" alt="Growth rate vs interaction index for CH4 and H2" width="600">

Both fuels become unstable above a threshold interaction index. The
thresholds and the subsequent behaviour differ markedly:

| Quantity | CH4/air | H2/air |
|---|---|---|
| tau [ms] | 26.64 | 4.99 |
| Instability threshold (sigma = 0) | n ~ 0.37 | n ~ 0.60 |
| sigma at n = 3 [1/s] | +58.4 | +167.5 |
| Frequency range over the sweep [Hz] | 86.9 - 97.1 | 87.6 - 114.6 |

Hydrogen crosses into instability at a *higher* interaction index than
methane in this configuration, but once unstable its growth rate rises
roughly three times faster. The curves cross near n ~ 0.85. Expressed in
engineering terms: hydrogen is not more prone to instability here, but the
instability it produces is considerably more violent.

<img src="../results/figures/fuel_comparison_frequency_vs_n.png" alt="Mode frequency along the sweep" width="600">

The frequency plot is the branch-tracking check: both curves start near
the acoustic-only eigenfrequency (85.75 Hz) and vary smoothly and
monotonically, confirming that the solver remained on the fundamental mode
throughout. Hydrogen shifts the mode frequency more strongly (up to
115 Hz) than methane (up to 97 Hz), consistent with its stronger coupling.

### 3.7 Rayleigh criterion verification

The Rayleigh index, proportional to Re(p'* x Q') at the flame location,
was computed at n = 2.0 for both calibrated flames and compared to the
sign of the growth rate obtained from the eigenvalue solver:

| Case | sigma [1/s] | Rayleigh index | Signs match |
|---|---|---|---|
| H2/air, n = 2.0 | +114.62 | +0.00066 | Yes |
| CH4/air, n = 2.0 | +46.65 | +0.00030 | Yes |

Both cases are unstable at n = 2, and in both the Rayleigh index is
positive: heat release is in phase with pressure, feeding the oscillation.
The eigenvalue solver and the Rayleigh criterion agree, as they must for a
linear system.

### 3.8 Limitations

Three assumptions limit the quantitative reach of these results:

- **No damping.** The open end is idealized (p = 0, no radiation
  impedance) and no viscous or thermal losses are modelled. With no
  dissipation to oppose it, the flame source term alone sets the stability
  boundary, which is why both fuels eventually become unstable. Realistic
  losses would raise both thresholds and could change their ordering.
- **Uniform mean temperature.** Both ducts use the same speed of sound
  and density, whereas in reality the temperature jumps from ~300 K to
  ~2230 K across the flame (a factor ~2.7 on c, ~7 on rho). This is a
  significant simplification of the acoustic field downstream of the flame.
- **n and d are not calibrated.** The interaction index is swept as a free
  parameter rather than derived from physics, and the convective length
  scale d entering tau is an order-of-magnitude estimate. The absolute
  position of the stability thresholds therefore carries no quantitative
  weight; only the comparison between the two fuels, computed with
  identical assumptions, is meaningful.


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
