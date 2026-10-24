"""
Eigenmode spectrum of a Rijke-type tube with a compact flame (closed-open
tube), for the CH4/air and H2/air flames of the YAML configurations.

The first acoustic modes are followed by continuation in the interaction
index n, each from its acoustic eigenfrequency f_j = (2j - 1) c / (4 L) at a
very small n. The flame gain n is the same at every frequency and the model
has no damping, so higher modes can become unstable before the fundamental.

Outputs: a table with the growth rate of each mode and its instability
threshold, the spectrum at a fixed n, and the growth rate of each mode
versus n.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

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

N_MODES = 6            
N_MIN = 0.005         
N_POINTS = 1500        
N_SPECTRUM = 1.0       

n_values = np.linspace(N_MIN, reference["sweep"]["n_max"], N_POINTS)

inlet = ClosedEnd()
outlet = OpenEnd()
duct_one = Duct(flame_position, speed_of_sound)
duct_two = Duct(tube_length - flame_position, speed_of_sound)



# Helpers


def acoustic_frequency(j):
    """Frequency [Hz] of acoustic mode j of the closed-open tube."""
    return (2 * j - 1) * speed_of_sound / (4 * tube_length)


def build_characteristic(tau):
    """characteristic(n) -> residual function of the tube with a flame of delay tau."""
    def characteristic(n):
        flame = CompactFlame(NTauFlame(n, tau), coupling_factor)
        network = ThermoAcousticNetwork(inlet, [duct_one, flame, duct_two], outlet)
        return network.characteristic_value
    return characteristic


def first_instability(n, sigma):
    """Smallest n at which sigma turns positive (None if it never does).
    Returns n[0] when the mode is already unstable at the smallest n."""
    if sigma[0] > 0:
        return n[0]
    upward = np.where(np.diff(np.sign(sigma)) > 0)[0]
    if not len(upward):
        return None
    i = upward[0]
    a = sigma[i] / (sigma[i] - sigma[i + 1])
    return n[i] + a * (n[i + 1] - n[i])



# Flame calibration and mode tracking


tau = {name: FlameCalibrator.from_config(cfg).calibrate(phi=cfg["fuel"]["equivalence_ratio"]).tau
       for name, cfg in configs.items()}

omega = {
    name: {j: track_mode(build_characteristic(tau[name]), n_values,
                         2 * np.pi * acoustic_frequency(j) + 1j)
           for j in range(1, N_MODES + 1)}
    for name in configs
}
sigma = {name: {j: StabilityAnalyzer.growth_rate(omega[name][j]) for j in omega[name]}
         for name in configs}

i_spectrum = int(np.argmin(abs(n_values - N_SPECTRUM)))
print(f"Interaction index of the spectrum: n = {n_values[i_spectrum]:.3f}")

for name in configs:
    print(f"\n{name}/air: tau = {tau[name] * 1e3:.2f} ms")
    print(f"  {'mode':>4} {'f0 [Hz]':>9} {'f(n=1) [Hz]':>12} {'sigma(n=1) [1/s]':>17}   first instability")
    thresholds = {}
    for j in range(1, N_MODES + 1):
        s = sigma[name][j]
        if not np.all(np.isfinite(s)):
            print(f"  {j:4d}  branch lost: the mode could not be followed over the whole sweep")
            continue
        n_c = first_instability(n_values, s)
        thresholds[j] = n_c
        if n_c is None:
            label = f"none up to n = {n_values[-1]:.1f}"
        elif s[0] > 0:
            label = f"unstable for all n >= {n_values[0]}"
        else:
            label = f"n = {n_c:.3f}"
        print(f"  {j:4d} {acoustic_frequency(j):9.2f} {omega[name][j][i_spectrum].real / (2 * np.pi):12.2f} "
              f"{s[i_spectrum]:17.1f}   {label}")
    unstable = {j: n_c for j, n_c in thresholds.items() if n_c is not None}
    if unstable:
        j_first = min(unstable, key=unstable.get)
        print(f"  First instability of the system: mode {j_first} "
              f"({acoustic_frequency(j_first):.0f} Hz) at n = {unstable[j_first]:.3f}")


# Figures

style = {"CH4": ("tab:blue", "o", "CH$_4$/air"), "H2": ("tab:orange", "D", "H$_2$/air")}

# Spectrum at a fixed interaction index
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axhline(0, color="black", ls="--", lw=1)
top = max(np.nanmax([sigma[name][j][i_spectrum] for j in sigma[name]]) for name in configs)
bottom = min(np.nanmin([sigma[name][j][i_spectrum] for j in sigma[name]]) for name in configs)
ax.axhspan(0, top * 1.25, color="red", alpha=0.07)
for name, (color, marker, label) in style.items():
    f_modes = [omega[name][j][i_spectrum].real / (2 * np.pi) for j in omega[name]]
    s_modes = [sigma[name][j][i_spectrum] for j in sigma[name]]
    ax.scatter(f_modes, s_modes, s=70, color=color, marker=marker, edgecolor="black",
               label=label, zorder=3)
    for j, (f, s) in enumerate(zip(f_modes, s_modes), start=1):
        ax.annotate(str(j), (f, s), textcoords="offset points", xytext=(7, 5), fontsize=8, color=color)
ax.text(0.01, 0.97, "Unstable ($\\sigma$ > 0)", transform=ax.transAxes, color="darkred", fontsize=9, va="top")
ax.set_ylim(bottom - 0.1 * abs(bottom) - 5, top * 1.25)
ax.set_xlabel("Frequency $f$ [Hz]")
ax.set_ylabel("Growth rate $\\sigma$ [1/s]")
ax.set_title(f"Eigenmode spectrum at $n$ = {n_values[i_spectrum]:.1f} (numbers = mode index)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(FIG_DIR / "mode_spectrum.png", dpi=150)

# Growth rate of each mode versus n
fig2, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), sharey=True)
colors = plt.cm.viridis(np.linspace(0.05, 0.9, N_MODES))
for ax, name in zip(axes, configs):
    for j in range(1, N_MODES + 1):
        ax.plot(n_values, sigma[name][j], color=colors[j - 1], lw=2.6 if j == 1 else 1.4,
                label=f"mode {j} ({acoustic_frequency(j):.0f} Hz)")
    ax.axhline(0, color="black", ls="--", lw=1)
    ax.set_title(f"{style[name][2]}  ($\\tau$ = {tau[name] * 1e3:.2f} ms)")
    ax.set_xlabel("Interaction index $n$ [-]")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("Growth rate $\\sigma$ [1/s]")
axes[1].legend(loc="upper left", fontsize=8)
fig2.tight_layout()
fig2.savefig(FIG_DIR / "mode_growth_vs_n.png", dpi=150)

plt.show()