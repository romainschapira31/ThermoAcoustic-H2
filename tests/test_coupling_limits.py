"""Analytical checks of the coupled solver.

For a closed-open tube of length L = L1 + L2 with a compact flame at x = L1,
the characteristic equation derived by hand from the transfer matrices is

    F(omega) = cos(kL) - xi * sin(k L1) * sin(k L2),
    xi = coupling_factor * n * exp(-i * omega * tau),   k = omega / c.

With xi = 0 it reduces to cos(kL) = 0, the acoustic eigenfrequencies. Two
exact consequences are used here as tests of the solver.
"""

import numpy as np
from scipy.optimize import newton

from thermoacoustic.acoustics import Duct, ClosedEnd, OpenEnd
from thermoacoustic.flame import NTauFlame, CompactFlame
from thermoacoustic.solver import ThermoAcousticNetwork

C = 343.0
L1, L2 = 0.25, 0.75
L = L1 + L2
COUPLING = 1.0
OMEGA0 = np.pi * C / (2 * L)       

TAU_H2 = 4.326e-3  
TAU_CH4 = 26.64e-3   


def build(n, tau, l1=L1):
    flame = CompactFlame(NTauFlame(n, tau), COUPLING)
    return ThermoAcousticNetwork(
        ClosedEnd(), [Duct(l1, C), flame, Duct(L - l1, C)], OpenEnd()
    )


def neutral_point(tau, m):
    """Exact neutral stability point.

    For real omega, cos(kL), sin(k L1) and sin(k L2) are real, so
    Im F = -Im(xi) sin(k L1) sin(k L2) = 0 forces sin(omega tau) = 0, i.e.
    omega = m pi / tau. Then xi = (-1)^m * coupling * n, and F = 0 gives n.
    """
    omega = m * np.pi / tau
    k = omega / C
    xi = np.cos(k * L) / (np.sin(k * L1) * np.sin(k * L2))
    return omega, (-1) ** m * xi / COUPLING


FUNDAMENTAL_NEUTRAL = [("H2", TAU_H2, 1), ("CH4", TAU_CH4, 5)]


def test_residual_vanishes_at_predicted_neutral_point():
    for _, tau, m in FUNDAMENTAL_NEUTRAL:
        omega_c, n_c = neutral_point(tau, m)
        assert n_c > 0
        assert abs(build(n_c, tau).characteristic_value(omega_c)) < 1e-9


def test_stability_changes_sign_at_predicted_threshold():
    """Slightly below n_c the mode decays, slightly above it grows."""
    for _, tau, m in FUNDAMENTAL_NEUTRAL:
        omega_c, n_c = neutral_point(tau, m)
        below = newton(build(0.98 * n_c, tau).characteristic_value, omega_c + 1j)
        above = newton(build(1.02 * n_c, tau).characteristic_value, omega_c + 1j)
        assert -below.imag < 0
        assert -above.imag > 0


def first_order_shift(n, tau, l1=L1):
    """First-order eigenfrequency shift of the fundamental for small n."""
    xi = COUPLING * n * np.exp(-1j * OMEGA0 * tau)
    k0 = OMEGA0 / C
    return -(C / L) * xi * np.sin(k0 * l1) * np.sin(k0 * (L - l1))


def relative_error(n, tau):
    omega = newton(build(n, tau).characteristic_value, OMEGA0 + 1j)
    delta_num = omega - OMEGA0
    delta_th = first_order_shift(n, tau)
    return abs(delta_num - delta_th) / abs(delta_th)


def test_first_order_shift_converges_linearly_in_n():
    """The error of a first-order expansion must scale as n (relative to the
    shift, which itself scales as n): dividing n by 10 divides it by ~10."""
    for tau in (TAU_H2, TAU_CH4):
        coarse = relative_error(0.01, tau)
        fine = relative_error(0.001, tau)
        assert fine < 1e-2
        assert 7 < coarse / fine < 13


def test_growth_rate_follows_flame_position():
    """At small n, sigma is proportional to sin(pi * L1 / L): the Rayleigh
    product p'u' of the quarter-wave mode peaks at L/2 in a closed-open tube."""
    n, tau = 1e-3, TAU_H2

    def sigma(l1):
        omega = newton(build(n, tau, l1).characteristic_value, OMEGA0 + 1j)
        return -omega.imag

    reference = sigma(L / 2)
    for l1 in (0.1, 0.25, 0.4):
        assert np.isclose(sigma(l1) / reference, np.sin(np.pi * l1 / L), rtol=1e-2)