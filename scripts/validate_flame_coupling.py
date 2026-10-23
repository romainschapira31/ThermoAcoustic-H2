"""
Thermoacoustic stability of a Rijke-type tube with a compact flame.

Computes the growth rate of the fundamental acoustic mode as a function of
the flame interaction index n, for CH4/air and H2/air flames whose time
delay tau is calibrated from Cantera laminar flame speeds.

Eigenfrequencies are obtained by direct root-finding on the characteristic
equation, with progressive continuation of the initial guess along n. The
real part of omega is tracked throughout to confirm that the solver stays
on the fundamental mode branch.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

from thermoacoustic.acoustics import Duct, OpenEnd, ClosedEnd
from thermoacoustic.flame import NTauFlame, CompactFlame, FlameCalibrator
from thermoacoustic.solver import ThermoAcousticNetwork


# Configuration

tube_length = 1.0            
flame_position = 0.25        
speed_of_sound = 343.0      
coupling_factor = 1.0
characteristic_distance = 0.01   

n_sweep = np.linspace(0.1, 3.0, 60)
omega_start = 2 * np.pi * 85.75 + 1j  

inlet = ClosedEnd()
outlet = OpenEnd()
duct_one = Duct(flame_position, speed_of_sound)
duct_two = Duct(tube_length - flame_position, speed_of_sound)


# Helpers


def sweep_growth_rate(tau_fixed, n_values, omega_guess=omega_start):
    """Compute sigma(n) and f(n) at fixed tau, by direct root-finding
    with continuation of the initial guess."""
    sigmas, freqs = [], []
    for n in n_values:
        network = ThermoAcousticNetwork(
            inlet,
            [duct_one, CompactFlame(NTauFlame(n, tau_fixed), coupling_factor), duct_two],
            outlet,
        )
        omega = newton(network.characteristic_value, omega_guess)
        omega_guess = omega
        sigmas.append(-omega.imag)
        freqs.append(omega.real / (2 * np.pi))
    return np.array(sigmas), np.array(freqs)


def rayleigh_index(omega, duct_before_flame, flame):
    """Rayleigh index Re(p'* x Q') at the flame location.
    Positive means heat release in phase with pressure (destabilizing)."""
    a_in, b_in = inlet.coefficients()
    if b_in != 0:
        p0, u0 = 1.0, -a_in / b_in
    else:
        p0, u0 = 0.0, 1.0
    state = duct_before_flame.transfer_matrix(omega) @ np.array([p0, u0], dtype=complex)
    p_flame, u_flame = state
    Q_prime = flame.flame_model.heat_release_response(omega) * u_flame
    return np.real(np.conj(p_flame) * Q_prime)


# Flame calibration from Cantera


calibrator_ch4 = FlameCalibrator(
    fuel='CH4', oxidizer='O2:1, N2:3.76', n=1.0,
    characteristic_distance=characteristic_distance, mechanism='gri30.yaml',
)
tau_ch4 = calibrator_ch4.calibrate(phi=1.0).tau

calibrator_h2 = FlameCalibrator(
    fuel='H2', oxidizer='O2:1, N2:3.76', n=1.0,
    characteristic_distance=characteristic_distance, mechanism='h2o2.yaml',
)
tau_h2 = calibrator_h2.calibrate(phi=0.9).tau

print(f"Calibrated time delays (d = {characteristic_distance*1000:.0f} mm):")
print(f"  CH4/air, phi=1.0 : tau = {tau_ch4*1000:6.2f} ms")
print(f"  H2/air,  phi=0.9 : tau = {tau_h2*1000:6.2f} ms")


# Growth rate vs interaction index


sigma_ch4, freq_ch4 = sweep_growth_rate(tau_ch4, n_sweep)
sigma_h2,  freq_h2  = sweep_growth_rate(tau_h2,  n_sweep)

print("\nGrowth rate over n in [0.1, 3.0]:")
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
    sigma = -omega.imag
    RI = rayleigh_index(omega, duct_one, flame)
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
fig.savefig("../results/figures/fuel_comparison_sigma_vs_n.png", dpi=150)

fig2, ax2 = plt.subplots()
ax2.plot(n_sweep, freq_ch4, color='tab:blue', lw=2, label='CH$_4$/air')
ax2.plot(n_sweep, freq_h2, color='tab:orange', lw=2, label='H$_2$/air')
ax2.axhline(85.75, color='grey', ls=':', lw=1.5, label='Acoustic-only mode')
ax2.set_xlabel('Interaction index $n$ [-]')
ax2.set_ylabel('Mode frequency $f$ [Hz]')
ax2.set_title('Mode frequency along the sweep (branch tracking check)')
ax2.legend()
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig("../results/figures/fuel_comparison_frequency_vs_n.png", dpi=150)

plt.show()