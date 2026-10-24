import numpy as np
import pytest
from scipy.optimize import newton

from thermoacoustic.acoustics import Duct, ClosedEnd, OpenEnd
from thermoacoustic.flame import NTauFlame, CompactFlame
from thermoacoustic.solver import ThermoAcousticNetwork
from thermoacoustic.stability import StabilityAnalyzer

OMEGA_START = 2 * np.pi * 85.75 + 1j


def build(n, tau):
    flame = CompactFlame(NTauFlame(n, tau), coupling_factor=1.0)
    network = ThermoAcousticNetwork(
        ClosedEnd(), [Duct(0.25, 343.0), flame, Duct(0.75, 343.0)], OpenEnd()
    )
    return network, flame


def test_rayleigh_sign_matches_growth_rate():
    """With no damping, the flame is the only energy source, so the sign
    of sigma must equal the sign of the Rayleigh index."""
    for n, tau in [(2.0, 0.00499), (2.0, 0.02664)]:
        network, flame = build(n, tau)
        omega = newton(network.characteristic_value, OMEGA_START)
        analyzer = StabilityAnalyzer(network)
        assert np.sign(analyzer.growth_rate(omega)) == np.sign(
            analyzer.rayleigh_index(omega, flame)
        )


def test_rayleigh_index_vanishes_without_flame_response():
    network, flame = build(0.0, 0.005)
    omega = newton(network.characteristic_value, OMEGA_START)
    assert abs(StabilityAnalyzer(network).rayleigh_index(omega, flame)) < 1e-12


def test_state_at_first_element_is_inlet_state():
    network, _ = build(1.0, 0.005)
    state = StabilityAnalyzer(network).state_at(500.0, network.elements[0])
    assert np.allclose(state, network.inlet_state())


def test_state_at_unknown_element_raises():
    network, _ = build(1.0, 0.005)
    with pytest.raises(ValueError):
        StabilityAnalyzer(network).state_at(500.0, Duct(1.0, 343.0))