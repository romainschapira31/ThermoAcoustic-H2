import numpy as np
from scipy.optimize import newton


class ThermoAcousticNetwork():
    """Assembles a sequence of acoustic elements between two boundary
       conditions and evaluates the characteristic equation of the system.

       A root of the characteristic value (as a function of omega) corresponds
       to an eigenfrequency of the network.
    """

    def __init__(self, inlet, elements, outlet):
        """Store the inlet boundary condition, the list of acoustic elements
        (traversed in order), and the outlet boundary condition."""

        self.inlet = inlet
        self.elements = elements
        self.outlet = outlet

    def inlet_state(self):
        """State [p, u] satisfying the inlet boundary condition."""

        a_in, b_in = self.inlet.coefficients()
        if b_in != 0:
            p, u = 1.0, -a_in / b_in
        else:
            p, u = 0.0, 1.0
        return np.array([p, u], dtype=complex)

    def characteristic_value(self, omega):
        """Compute the residual of the outlet boundary condition for a state
        that already satisfies the inlet condition, propagated through all
        elements at the given angular frequency omega.

        A residual close to zero indicates that omega is close to an
        eigenfrequency of the network.
        """

        a_out, b_out = self.outlet.coefficients()

        state = self.inlet_state()

        for element in self.elements:
            state = element.transfer_matrix(omega) @ state

        p_out, u_out = state

        return a_out * p_out + b_out * u_out


def track_mode(characteristic, parameters, omega_start, max_step=10.0):
    """Follow one eigenmode along a parameter sweep.

    characteristic(p) returns the function omega -> residual of the network at
    parameter value p. The root at the first value is searched from
    omega_start; each following root is searched from the previous one.

    A root that moves by more than max_step [rad/s] between two consecutive
    parameter values is taken to be a jump onto another branch: the mode is
    then no longer followed and the remaining entries are NaN, as they are
    when Newton's method fails to converge.
    """
    omegas = np.full(len(parameters), complex(np.nan, np.nan))
    guess = omega_start
    for i, p in enumerate(parameters):
        try:
            root = newton(characteristic(p), guess, maxiter=100)
        except RuntimeError:
            break
        if i > 0 and abs(root - guess) > max_step:
            break
        omegas[i] = root
        guess = root
    return omegas