"""Stability analysis of a thermoacoustic network."""

import numpy as np


class StabilityAnalyzer:
    """Growth rate and Rayleigh index of an eigenmode of a network."""

    def __init__(self, network):
        """Store the network to analyse."""

        self.network = network

    @staticmethod
    def growth_rate(omega):
        """Growth rate sigma = -Im(omega); sigma > 0 means unstable."""

        return -omega.imag

    def is_unstable(self, omega):
        """Return True if the mode with eigenfrequency omega grows in time."""

        return self.growth_rate(omega) > 0

    def state_at(self, omega, element):
        """Acoustic state [p, u] at the inlet of the given element."""

        state = self.network.inlet_state()
        for current in self.network.elements:
            if current is element:
                return state
            state = current.transfer_matrix(omega) @ state
        raise ValueError("Element is not part of the network")

    def rayleigh_index(self, omega, flame):
        """Rayleigh index Re(p'* x Q') at the flame location.
        Positive means heat release in phase with pressure (destabilizing)."""

        p, u = self.state_at(omega, flame)
        q_prime = flame.flame_model.heat_release_response(omega) * u
        return np.real(np.conj(p) * q_prime)