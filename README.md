# ThermoAcoustic-H2

Low-order thermoacoustic stability modeling for premixed CH4/air and H2/air flames
in a Rijke-tube-type configuration. The network acoustic model (n-tau flame response)
is calibrated against Cantera laminar flame calculations to compare the stability
behavior of methane and hydrogen combustion.

## Context

This project was developed as a personal research project to explore the coupling
between acoustics and combustion, with a focus on the stability implications of
hydrogen as an alternative fuel for aviation gas turbines.

## Architecture

- `src/thermoacoustic/acoustics.py` — acoustic network elements (ducts, junctions, boundaries)
- `src/thermoacoustic/flame.py` — flame response models (n-tau, Cantera-based calibration)
- `src/thermoacoustic/solver.py` — network assembly and eigenvalue solver
- `src/thermoacoustic/stability.py` — Rayleigh criterion and fuel comparison tools

## Installation

```bash
pip install -e .
```

## Usage

```bash
python scripts/validate_acoustics.py
python main.py --config configs/ch4_rijke.yaml
```

## Results

(to be completed as the project progresses)

## Limitations

This is a low-order model, not a high-fidelity LES — see `docs/formulation.md` for
the full list of assumptions.
