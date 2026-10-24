# ThermoAcoustic-H2

Low-order thermoacoustic stability modeling of a Rijke-type tube with a
compact premixed flame, used to compare CH4/air and H2/air. The flame
response is described by Crocco's n-tau model, with the time delay tau
calibrated from Cantera laminar flame speeds.

## Context

This is a personal research project exploring the coupling between acoustics
and combustion, with a focus on what a switch to hydrogen could mean for
combustion instabilities in aviation gas turbines. It is a low-order study,
not a high-fidelity simulation: its value lies in the physical formulation,
the verification and validation of each building block, and the quality of
the code.

## Architecture

- `src/thermoacoustic/acoustics.py` — acoustic network elements (ducts, boundary conditions)
- `src/thermoacoustic/flame.py` — n-tau flame model, compact flame element, and `FlameCalibrator` (tau from Cantera)
- `src/thermoacoustic/solver.py` — network assembly and characteristic equation
- `src/thermoacoustic/stability.py` — growth rate and Rayleigh index analysis
- `src/thermoacoustic/config.py` — YAML configuration loading and validation
- `configs/` — one YAML file per fuel case
- `scripts/` — validation and study scripts
- `tests/` — unit tests and analytical checks
- `docs/formulation.md` — physical formulation, verification details, limitations

## Installation

```bash
pip install -e ".[dev]"
```

## Usage

```bash
# Flame calibration (Cantera) — validates flame speed and adiabatic temperature
python scripts/cantera_reference_case.py

# Acoustic-only verification — eigenfrequencies of a closed-open duct
python scripts/validate_acoustics.py

# Thermoacoustic coupling — CH4 vs H2 growth rate comparison and Rayleigh check
python scripts/validate_flame_coupling.py

# Sensitivity of the fuel comparison to the convective length d
python scripts/phase_sensitivity.py
# Sensitivity of the computed flame speeds to numerical and transport settings
python scripts/transport_check.py
```

The two study scripts (`validate_flame_coupling.py` and
`phase_sensitivity.py`) read their parameters (fuel, geometry, flame model,
sweep range) from the YAML files in `configs/`. The two reference scripts
reproduce fixed cases (analytical duct modes, published flame speed) and
keep their parameters in the script by design.

## Testing

```bash
pytest tests/ -v
```

The test suite covers the acoustic elements and boundary conditions, the
solver (first eigenfrequency against the analytical closed-open duct value),
the YAML configuration loader and its consistency checks, the stability
analysis (the sign of the Rayleigh index must match the sign of the growth
rate), and analytical checks of the coupled solver: the exact stability
threshold, the first-order frequency shift, and the dependence on flame
position.

## Results

### Flame calibration (Cantera)

Laminar flame speed and adiabatic flame temperature of a stoichiometric
CH4/air mixture, computed with GRI-Mech 3.0 at T = 298.15 K and
p = 101325 Pa, matching the conditions of the literature references.

| Quantity | This work | Reference |
|---|---|---|
| Laminar flame speed, CH4/air | 0.375 m/s | 0.36 ± 0.017 m/s (compilation of twelve studies, FFCM-1) |
| Laminar flame speed, H2/air | 2.311 m/s | 2.11 ± 0.11 m/s (compilation of nine studies, FFCM-2) |
| Adiabatic flame temperature, CH4/air | 2229.5 K | 2224.3 K (same mechanism and software, arXiv:2503.11826) |

> The reference flame speeds are consensus values compiled by the Stanford
> FFCM initiative from twelve (CH4) and nine (H2) experimental studies. The
> original papers were not accessed; see `docs/formulation.md`, References
> 1 and 5.

The computed methane flame speed is 4.3% above the compiled measurement
(+1.5 cm/s), within its uncertainty. The hydrogen flame speed is 9.5% above
(+20 cm/s), outside it. The gap is not numerical: mesh refinement and domain
size change it by less than 0.5%, and the full GRI-Mech 3.0 gives the same
value. It comes from the transport model: with multicomponent transport and
thermal diffusion, the computed value is 1.1% below the measurement
(`scripts/transport_check.py`). The project keeps Cantera's default
transport for both fuels. The adiabatic temperature differs by 0.23% from an
independent calculation using the same mechanism and software, which checks
the implementation rather than the mechanism. Grid independence was verified
through a convergence study (34 → 139 → 253 points); see
`docs/formulation.md` for the full table.

The flame-speed uncertainty moves the stability thresholds by tens of
percent (`docs/formulation.md`, Section 3.8), but at d = 10 mm hydrogen
still crosses later than methane for all the flame speeds tested.

### Acoustic network verification

The transfer-matrix acoustic solver was verified against the analytical
eigenfrequencies of a closed-open duct (L = 1 m, c = 343 m/s): all 6
computed modes agree with theory within 0.3%, the residual being the 1 Hz
resolution of the frequency sweep.

### Thermoacoustic coupling and fuel comparison

A compact flame element (Crocco's n-tau model) is coupled to the acoustic
network, making the eigenvalue problem complex: each mode carries a growth
rate sigma = -Im(omega), with sigma > 0 indicating instability.

The time delay is estimated as a convective transit time, tau = d / S_L,
with S_L from Cantera and d = 10 mm, for both fuels at phi = 1.0:

| Fuel | S_L [m/s] | tau [ms] |
|---|---|---|
| CH4/air | 0.375 | 26.64 |
| H2/air | 2.311 | 4.33 |

The ratio of the two delays (6.2) is fixed by the chemistry. The
interaction index n and the length d are **not** calibrated.

**Exact stability threshold.** A neutral mode requires omega·tau = m·pi,
which gives the critical interaction index in closed form. The solver
reproduces it to the digits shown:

| Fuel | f at threshold [Hz] | n_c (exact) | n_c (solver) |
|---|---|---|---|
| CH4/air | 93.8 | 0.369 | 0.369 |
| H2/air | 115.6 | 1.029 | 1.029 |

Instability sets in when the flame has pulled the mode frequency up to
omega·tau = m·pi, starting from the acoustic fundamental at 85.75 Hz.

Sweeping n at each fuel's calibrated tau gives:

<img src="results/figures/fuel_comparison_sigma_vs_n.png" alt="Growth rate vs interaction index for CH4 and H2" width="600">

Both fuels become unstable above their threshold. At d = 10 mm, hydrogen
crosses later (n = 1.03 against 0.37 for methane) but its growth rate then
rises faster: it overtakes methane's above n ≈ 1.3 and reaches +169 against
+58 1/s at n = 3. The frequency of the tracked mode stays smooth and
monotonic along both sweeps, confirming that the solver remains on the
fundamental mode:

<img src="results/figures/fuel_comparison_frequency_vs_n.png" alt="Mode frequency along the sweep" width="600">

**What does not carry over.** These statements hold for d = 10 mm. The
order in which the two fuels become unstable depends on d: at d = 7.5 mm
methane is unstable at small n while hydrogen is not, and at d = 15 mm the
reverse holds. With d uncalibrated, this model does not support a ranking of
the two fuels. What is robust is the delay ratio set by the chemistry, the
phase mechanism, and the exact threshold formula. A ranking would require a
calibrated d or a measured flame transfer function (see
`docs/formulation.md`, Section 3.8).

The Rayleigh index agrees in sign with the growth rate for both fuels, as it
must when the flame is the only source of acoustic energy.

**Methodological note.** An initial attempt at a 2D (n, tau) stability map
was discarded. Newton's method on this transcendental characteristic
equation does not guarantee branch tracking, and the map silently drifted
onto a spurious mode. The error was found by cross-checking against direct
1D sweeps, which are used for all results above. See
`docs/formulation.md` (Section 3.4) for the full diagnosis.

## Limitations

This is a low-order model, and its quantitative reach is limited:

- **No damping.** The open end is idealized and no losses are modelled, so
  the flame alone sets the stability boundary. Realistic losses would raise
  both thresholds and could change their ordering.
- **Uniform mean temperature.** Both ducts share the same speed of sound,
  whereas the temperature rises by a factor of ~7.5 across a real flame.
- **Uncalibrated parameters.** The interaction index n, the coupling factor
  (only their product matters) and the convective length d are not derived
  from physics. Absolute thresholds carry no quantitative weight, and the
  fuel comparison depends on d.
- **Flame speed from a simple transport model, one equivalence ratio.** With
  Cantera's default transport, the hydrogen flame speed is 9.5% above the
  measurement; a more complete transport model removes the gap. Only
  phi = 1.0 is reported.
- **No external validation of the coupled model.** It is verified against
  analytical limits but not compared with a published or measured
  thermoacoustic case.

See `docs/formulation.md` for the derivation, the calibration, and the
complete list of assumptions.