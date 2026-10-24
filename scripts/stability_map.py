"""
Stability map of the fundamental mode in the (n, tau) plane, for a Rijke-type
tube (closed-open) with a compact flame.

Each row of the map is a 1D continuation along n at fixed tau, started from
the acoustic fundamental at a very small n, where the mode is known to sit
next to its acoustic value. Rows are independent: no root is carried from one
tau to the next, which is what made an earlier 2D map jump onto spurious
branches (see docs/formulation.md, Section 3.4).

The zero contour of the map is compared with the exact neutral points
omega * tau = m * pi (Section 3.6), and the CH4/air and H2/air flames of the
YAML configurations are placed on the map, with the range of tau covered when
the convective length d varies between half and twice its configured value.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from thermoacoustic.acoustics import Duct, OpenEnd, ClosedEnd
from thermoacoustic.config import load_config, check_same_setup
from thermoacoustic.flame import NTauFlame, CompactFlame, FlameCalibrator
from thermoacoustic.solver import ThermoAcousticNetwork, track_mode
from thermoacoustic.stability import StabilityAnalyzer


# Configuration (read from YAML)


ROOT = Path(__file__).resolve().parent.parent   
CONFIG_DIR = ROOT / "configs"
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

configs = {
    "CH4": load_config(CONFIG_DIR / "ch4_rijke.yaml"),
    "H2": load_config(CONFIG_DIR / "h2_rijke.yaml"),
}
check_same_setup(*configs.values())   

reference = configs["CH4"]
tube_length = reference["geometry"]["tube_length"]
flame_position = reference["geometry"]["flame_position"]
speed_of_sound = reference["geometry"]["speed_of_sound"]
coupling_factor = reference["flame_model"]["coupling_factor"]
d_reference = reference["flame_model"]["characteristic_distance"]
n_max = reference["sweep"]["n_max"]

N_MIN = 0.005                   
N_POINTS = 300                  
TAU_MIN, TAU_MAX = 1e-3, 60e-3  
N_TAU = 200                    
D_FACTORS = (0.5, 2.0)       

n_values = np.linspace(N_MIN, n_max, N_POINTS)
tau_values = np.linspace(TAU_MIN, TAU_MAX, N_TAU)

inlet = ClosedEnd()
outlet = OpenEnd()
duct_one = Duct(flame_position, speed_of_sound)
duct_two = Duct(tube_length - flame_position, speed_of_sound)

omega_acoustic = np.pi * speed_of_sound / (2 * tube_length)   


# Helpers


def build_characteristic(tau):
    """characteristic(n) -> residual function of the tube with a flame of delay tau."""
    def characteristic(n):
        flame = CompactFlame(NTauFlame(n, tau), coupling_factor)
        network = ThermoAcousticNetwork(inlet, [duct_one, flame, duct_two], outlet)
        return network.characteristic_value
    return characteristic


def exact_threshold(tau, m):
    """Critical interaction index of the neutral mode omega * tau = m * pi."""
    k = m * np.pi / (speed_of_sound * tau)
    return ((-1) ** m * np.cos(k * tube_length)
            / (coupling_factor * np.sin(k * flame_position) * np.sin(k * (tube_length - flame_position))))





REFINEMENTS = (1, 2, 4, 8)


def follow_row(tau):
    """Follow the fundamental along n at fixed tau, on the grid n_values.

    If the branch is lost (a step is too coarse for the root to be followed),
    the row is repeated with finer steps and sampled back onto n_values.
    Returns the roots and the refinement that was needed."""
    for refinement in REFINEMENTS:
        n_fine = np.linspace(N_MIN, n_max, (N_POINTS - 1) * refinement + 1)
        roots = track_mode(build_characteristic(tau), n_fine, omega_acoustic + 1j)
        if np.all(np.isfinite(roots)):
            break
    return roots[::refinement], refinement


rows = [follow_row(tau) for tau in tau_values]
omega = np.array([roots for roots, _ in rows])
refined = sum(refinement > 1 for _, refinement in rows)
sigma = StabilityAnalyzer.growth_rate(omega)
frequency = omega.real / (2 * np.pi)
tracked = np.all(np.isfinite(sigma), axis=1)

print(f"Rows followed over the whole sweep: {tracked.sum()} / {N_TAU} ({refined} needed finer steps)")
print(f"Frequency of the tracked mode: {np.nanmin(frequency):.1f} - {np.nanmax(frequency):.1f} Hz")


# Comparison of the zero contour with the exact neutral points

errors, exact_points = [], []
for i, tau in enumerate(tau_values):
    s = sigma[i]
    if not tracked[i] or s[0] > 0:
        continue                                   
    upward = np.where(np.diff(np.sign(s)) > 0)[0]
    if not len(upward):
        continue
    j = upward[0]
    a = s[j] / (s[j] - s[j + 1])
    n_numerical = n_values[j] + a * (n_values[j + 1] - n_values[j])
    omega_c = omega[i, j].real + a * (omega[i, j + 1].real - omega[i, j].real)
    m = int(round(omega_c * tau / np.pi))
    n_exact = exact_threshold(tau, m)
    errors.append(abs(n_numerical - n_exact))
    exact_points.append((n_exact, tau))

print(f"Zero crossings compared with the exact formula: {len(errors)}")
if errors:
    print(f"  |n_numerical - n_exact|: median = {np.median(errors):.5f}, "
          f"95th percentile = {np.percentile(errors, 95):.5f}, max = {max(errors):.4f}")
exact_points = np.array(exact_points)


# Flame calibration and figure


tau_fuel = {name: FlameCalibrator.from_config(cfg).calibrate(phi=cfg["fuel"]["equivalence_ratio"]).tau
            for name, cfg in configs.items()}
style = {"CH4": ("tab:blue", "CH$_4$/air"), "H2": ("tab:orange", "H$_2$/air")}
d_mm = d_reference * 1e3

fig, ax = plt.subplots(figsize=(8.2, 6.0))
vmax = 100
levels = np.linspace(-vmax, vmax, 21)
cs = ax.contourf(n_values, tau_values * 1e3, np.ma.masked_invalid(np.clip(sigma, -vmax, vmax)),
                 levels=levels, cmap="RdBu_r")
ax.contour(n_values, tau_values * 1e3, np.ma.masked_invalid(sigma), levels=[0],
           colors="black", linewidths=1.6)
if len(exact_points):
    shown = exact_points[::4]
    ax.plot(shown[:, 0], shown[:, 1] * 1e3, "o", ms=3.6, mfc="gold", mec="black", mew=0.6, ls="none")

handles = [Line2D([], [], color="black", lw=1.6, label="Solver: $\\sigma$ = 0"),
           Line2D([], [], marker="o", ls="none", mfc="gold", mec="black", ms=5,
                  label="Exact threshold ($\\omega\\tau = m\\pi$)")]
for name, (color, label) in style.items():
    t = tau_fuel[name] * 1e3
    ax.axhline(t, color=color, lw=2.2)
    ax.fill_betweenx([D_FACTORS[0] * t, D_FACTORS[1] * t], n_max - 0.09, n_max, color=color, alpha=0.95)
    handles.append(Line2D([], [], color=color, lw=2.2,
                          label=f"{label}, d = {d_mm:.0f} mm (bar: d = {D_FACTORS[0] * d_mm:.0f}-{D_FACTORS[1] * d_mm:.0f} mm)"))

ax.set_xlim(0, n_max)
ax.set_ylim(TAU_MIN * 1e3, TAU_MAX * 1e3)
ax.set_xlabel("Interaction index $n$ [-]")
ax.set_ylabel("Time delay $\\tau$ [ms]")
ax.set_title("Fundamental mode: growth rate in the ($n$, $\\tau$) plane")
fig.colorbar(cs, ax=ax, label="Growth rate $\\sigma$ [1/s] (clipped at $\\pm$100)")
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig(FIG_DIR / "stability_map.png", dpi=150)

plt.show()