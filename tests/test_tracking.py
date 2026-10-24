import numpy as np

from thermoacoustic.acoustics import Duct, ClosedEnd, OpenEnd
from thermoacoustic.flame import NTauFlame, CompactFlame
from thermoacoustic.solver import ThermoAcousticNetwork, track_mode

C = 343.0
L1, L2 = 0.25, 0.75
COUPLING = 1.0
OMEGA0 = np.pi * C / (2 * (L1 + L2))
TAU_H2 = 4.326e-3      # s, hydrogen delay at phi = 1.0


def characteristic(tau):
    """characteristic(n) -> residual function of the closed-open tube with a flame."""
    def at(n):
        flame = CompactFlame(NTauFlame(n, tau), COUPLING)
        network = ThermoAcousticNetwork(ClosedEnd(), [Duct(L1, C), flame, Duct(L2, C)], OpenEnd())
        return network.characteristic_value
    return at


def test_tracking_reproduces_the_exact_threshold():
    """Following the fundamental from n ~ 0 must cross sigma = 0 at the n_c given
    by the closed-form neutral point omega * tau = pi."""
    n = np.linspace(0.005, 1.5, 600)
    omega = track_mode(characteristic(TAU_H2), n, OMEGA0 + 1j)
    sigma = -omega.imag
    i = np.where(np.diff(np.sign(sigma)) > 0)[0][0]
    a = sigma[i] / (sigma[i] - sigma[i + 1])
    n_numerical = n[i] + a * (n[i + 1] - n[i])

    k = np.pi / (C * TAU_H2)
    n_exact = -np.cos(k * (L1 + L2)) / (np.sin(k * L1) * np.sin(k * L2))
    assert abs(n_numerical - n_exact) < 2e-3


def test_tracking_stops_at_a_jump():
    """A root that jumps discontinuously must not be followed across the jump."""
    def jumping(p):
        target = 0.0 if p < 0.5 else 100.0
        return lambda omega: omega - target

    p = np.linspace(0.0, 1.0, 11)
    omega = track_mode(jumping, p, 0.0 + 0j, max_step=10.0)
    assert np.all(np.isfinite(omega[p < 0.5]))
    assert np.all(np.isnan(omega[p >= 0.5]))