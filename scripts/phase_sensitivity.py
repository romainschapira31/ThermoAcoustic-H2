"""
Small-n stability of the acoustic fundamental versus the convective length d.

To first order in the interaction index n, the growth rate of the
fundamental of the closed-open tube is

    sigma ~ -(c / 2L) * sin(pi * L1 / L) * cf * n * sin(omega0 * tau),
    tau = d / S_L,   omega0 = pi * c / (2 L),

so the mode is unstable for arbitrarily small n whenever sin(omega0 * tau) < 0.
Since tau depends on d, so does the stability of each fuel. This script
tabulates omega0 * tau (mod 2 pi) and its sine for both fuels over a range of
d, using the geometry and fuels of the YAML configurations.
"""

from pathlib import Path

import numpy as np

from thermoacoustic.config import load_config, check_same_setup
from thermoacoustic.flame import FlameCalibrator

ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = ROOT / "configs"

configs = {
    "CH4": load_config(CONFIG_DIR / "ch4_rijke.yaml"),
    "H2": load_config(CONFIG_DIR / "h2_rijke.yaml"),
}
check_same_setup(*configs.values())

geometry = configs["CH4"]["geometry"]
omega0 = np.pi * geometry["speed_of_sound"] / (2 * geometry["tube_length"])


def flame_speed(cfg):
    """Laminar flame speed S_L [m/s] of a fuel configuration, from Cantera."""
    fuel = cfg["fuel"]
    calibrator = FlameCalibrator(
        fuel=fuel["name"],
        oxidizer=fuel["oxidizer"],
        characteristic_distance=cfg["flame_model"]["characteristic_distance"],
        mechanism=fuel["mechanism"],
        T_in=fuel["T_fresh"],
        P_in=fuel["pressure"],
    )
    return calibrator.compute_flame_velocity(fuel["equivalence_ratio"])


S_L = {name: flame_speed(cfg) for name, cfg in configs.items()}

print(f"omega0 = {omega0:.2f} rad/s  (f0 = {omega0 / (2 * np.pi):.2f} Hz)")
for name, cfg in configs.items():
    print(f"{name}: phi = {cfg['fuel']['equivalence_ratio']}, S_L = {S_L[name]:.4f} m/s")

print(f"\n{'d [mm]':>7} {'fuel':>5} {'tau [ms]':>9} {'w0*tau mod 2pi':>15} {'sin':>7}  small-n regime")
for d in (0.005, 0.0075, 0.010, 0.015, 0.020):
    for name in configs:
        tau = d / S_L[name]
        phase = (omega0 * tau) % (2 * np.pi)
        regime = "stable" if np.sin(phase) > 0 else "unstable"
        print(f"{d * 1e3:7.1f} {name:>5} {tau * 1e3:9.2f} {phase:15.2f} {np.sin(phase):7.2f}  {regime}")