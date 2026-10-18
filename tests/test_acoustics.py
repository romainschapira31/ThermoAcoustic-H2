import numpy as np
import pytest

from thermoacoustic.acoustics import Duct, ClosedEnd, OpenEnd
from thermoacoustic.solver import ThermoAcousticNetwork


# Duct
def test_duct_transfer_matrix_zero_frequency():
    """At omega=0, the transfer matrix should reduce to the identity matrix
    (no propagation effect at zero frequency)."""
    duct = Duct(length=1.0, speed_of_sound=343.0)
    T = duct.transfer_matrix(omega=0.0)
    assert np.allclose(T, np.eye(2), atol=1e-8)


def test_duct_transfer_matrix_shape_and_dtype():
    """The transfer matrix must be a 2x2 complex array."""
    duct = Duct(length=1.0, speed_of_sound=343.0)
    T = duct.transfer_matrix(omega=100.0)
    assert T.shape == (2, 2)
    assert np.iscomplexobj(T)


# Boundary conditions
def test_closed_end_coefficients():
    a, b = ClosedEnd().coefficients()
    assert (a, b) == (0.0, 1.0)


def test_open_end_coefficients():
    a, b = OpenEnd().coefficients()
    assert (a, b) == (1.0, 0.0)


# ThermoAcousticNetwork 
def test_first_eigenfrequency_closed_open_tube():
    """The first eigenfrequency of a closed-open duct should match the
    analytical formula f_n = (2n-1)*c/(4L), within 1% for n=1."""
    L = 1.0
    c = 343.0

    network = ThermoAcousticNetwork(
        inlet=ClosedEnd(),
        elements=[Duct(length=L, speed_of_sound=c)],
        outlet=OpenEnd(),
    )

    f_analytic = c / (4 * L)  # first mode, n=1

    frequencies = np.linspace(f_analytic - 5, f_analytic + 5, 2000)
    residues = [
        abs(network.characteristic_value(2 * np.pi * f)) for f in frequencies
    ]
    f_found = frequencies[int(np.argmin(residues))]

    relative_error = abs(f_found - f_analytic) / f_analytic
    assert relative_error < 0.01


def test_characteristic_value_is_complex():
    """characteristic_value should return a complex number (even if its
    imaginary part is often negligible for a purely acoustic network)."""
    network = ThermoAcousticNetwork(
        inlet=ClosedEnd(),
        elements=[Duct(length=1.0, speed_of_sound=343.0)],
        outlet=OpenEnd(),
    )
    value = network.characteristic_value(omega=100.0)
    assert isinstance(value, complex) or isinstance(value, np.complexfloating)