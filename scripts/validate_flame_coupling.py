"""
Thermoacoustic stability of a Rijke-type tube with a compact flame.

Computes the growth rate of the fundamental acoustic mode as a function of
the flame interaction index n, for CH4/air and H2/air flames whose time
delay tau is calibrated from Cantera laminar flame speeds.

The fundamental mode is followed by continuation along n, started from the
acoustic eigenfrequency of the tube; if the root jumps onto another branch,
the script stops with an error. The frequency of the tracked mode is plotted
as a branch-tracking check, and the Rayleigh index is compared with the sign
of the growth rate.

All physical and numerical parameters are read from the YAML files in
configs/.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

from thermoacoustic.acoustics import Duct, OpenEnd, ClosedEnd
from thermoacoustic.config import load_config, check_same_setup
from thermoacoustic.flame import NTauFlame, CompactFlame, FlameCalibrator
from thermoacoustic.solver import ThermoAcousticNetwork, track_mode
from thermoacoustic.stability import StabilityAnalyzer


# Configuration (read from YAML)

ROOT = Path(__file__).resolve().parent.parent   # scripts/ -> project root
CONFIG_DIR = ROOT / "configs"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

cfg_ch4 = load_config(CONFIG_DIR / "ch4_rijke.yaml")
cfg_h2 = load_config(CONFIG_DIR / "h2_rijke.yaml")
check_same_setup(cfg_ch4, cfg_h2)   # geometry, numerics and phi must match

tube_length = cfg_ch4["geometry"]["tube_length"]
flame_position = cfg_ch4["geometry"]["flame_position"]
speed_of_sound = cfg_ch4["geometry"]["speed_of_sound"]
coupling_factor = cfg_ch4["flame_model"]["coupling_factor"]
characteristic_distance = cfg_ch4["flame_model"]["characteristic_distance"]

sweep = cfg_ch4["sweep"]
n_sweep = np.linspace(sweep["n_min"], sweep["n_max"], sweep["n_points"])

# Acoustic fundamental of the closed-open tube, where the root search starts
omega_acoustic = np.pi * speed_of_sound / (2 * tube_length)
omega_start = omega_acoustic + 1j

inlet = ClosedEnd()
outlet = OpenEnd()
duct_one = Duct(flame_position, speed_of_sound)
duct_two = Duct(tube_length - flame_position, speed_of_sound)


# Helpers

def sweep_growth_rate(tau_fixed, n_values):
    """Growth rate sigma(n) and frequency f(n) of the fundamental at fixed tau,
    by continuation along n. Raises if the branch is lost."""
    def characteristic(n):
        flame = CompactFlame(NTauFlame(n, tau_fixed), coupling_factor)
        network = ThermoAcousticNetwork(inlet, [duct_one, flame, duct_two], outlet)
        return network.characteristic_value

    omega = track_mode(characteristic, n_values, omega_start)
    if not np.all(np.isfinite(omega)):
        raise RuntimeError("The fundamental mode could not be followed over the whole sweep")
    return StabilityAnalyzer.growth_rate(omega), omega.real / (2 * np.pi)


def calibrate_tau(cfg):
    """Compute the flame time delay tau = d / S_L from a fuel configuration."""
    calibrator = FlameCalibrator.from_config(cfg)
    return calibrator.calibrate(phi=cfg["fuel"]["equivalence_ratio"]).tau


# Flame calibration from Cantera

tau_ch4 = calibrate_tau(cfg_ch4)
tau_h2 = calibrate_tau(cfg_h2)

print(f"Calibrated time delays (d = {characteristic_distance*1000:.0f} mm):")
print(f"  CH4/air, phi={cfg_ch4['fuel']['equivalence_ratio']} : tau = {tau_ch4*1000:6.2f} ms")
print(f"  H2/air,  phi={cfg_h2['fuel']['equivalence_ratio']} : tau = {tau_h2*1000:6.2f} ms")


# Growth rate vs interaction index

sigma_ch4, freq_ch4 = sweep_growth_rate(tau_ch4, n_sweep)
sigma_h2,  freq_h2  = sweep_growth_rate(tau_h2,  n_sweep)

print(f"\nGrowth rate over n in [{sweep['n_min']}, {sweep['n_max']}]:")
print(f"  CH4: sigma {sigma_ch4.min():8.1f} to {sigma_ch4.max():8.1f} 1/s"
      f"   | f {freq_ch4.min():.1f}-{freq_ch4.max():.1f} Hz")
print(f"  H2 : sigma {sigma_h2.min():8.1f} to {sigma_h2.max():8.1f} 1/s"
      f"   | f {freq_h2.min():.1f}-{freq_h2.max():.1f} Hz")


# Rayleigh criterion check

for label, tau_val, n_val in [("H2 ", tau_h2, 2.0),
                              ("CH4", tau_ch4, 2.0)]:
    flame = CompactFlame(NTauFlame(n_val, tau_val), coupling_factor)
    network = ThermoAcousticNetwork(inlet, [duct_one, flame, duct_two], outlet)
    omega = newton(network.characteristic_value, omega_start)
    analyzer = StabilityAnalyzer(network)
    sigma = analyzer.growth_rate(omega)
    RI = analyzer.rayleigh_index(omega, flame)
    print(f"{label}: n={n_val}, sigma={sigma:8.2f}, Rayleigh index={RI:+.5f}, "
          f"signs match: {np.sign(sigma) == np.sign(RI)}")


# Figures

fig, ax = plt.subplots()
ax.plot(n_sweep, sigma_ch4, color='tab:blue', lw=2,
        label=f'CH$_4$/air  ($\\tau$ = {tau_ch4*1000:.1f} ms)')
ax.plot(n_sweep, sigma_h2, color='tab:orange', lw=2,
        label=f'H$_2$/air  ($\\tau$ = {tau_h2*1000:.1f} ms)')
ax.axhline(0, color='black', ls='--', lw=1)
ax.fill_between(n_sweep, 0, max(sigma_h2.max(), sigma_ch4.max()) * 1.1,
                color='red', alpha=0.07)
ax.text(0.15, max(sigma_h2.max(), sigma_ch4.max()) * 0.9,
        'Unstable ($\\sigma$ > 0)', color='darkred', fontsize=9)
ax.set_xlabel('Interaction index $n$ [-]')
ax.set_ylabel('Growth rate $\\sigma$ [1/s]')
ax.set_title('Fundamental mode growth rate — CH$_4$ vs H$_2$')
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG_DIR / "fuel_comparison_sigma_vs_n.png", dpi=150)

fig2, ax2 = plt.subplots()
ax2.plot(n_sweep, freq_ch4, color='tab:blue', lw=2, label='CH$_4$/air')
ax2.plot(n_sweep, freq_h2, color='tab:orange', lw=2, label='H$_2$/air')
ax2.axhline(omega_acoustic / (2 * np.pi), color='grey', ls=':', lw=1.5, label='Acoustic-only mode')
ax2.set_xlabel('Interaction index $n$ [-]')
ax2.set_ylabel('Mode frequency $f$ [Hz]')
ax2.set_title('Mode frequency along the sweep (branch tracking check)')
ax2.legend()
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig(FIG_DIR / "fuel_comparison_frequency_vs_n.png", dpi=150)

plt.show()