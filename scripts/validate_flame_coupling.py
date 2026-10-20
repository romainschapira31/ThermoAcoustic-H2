import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from thermoacoustic.acoustics import Duct, OpenEnd, ClosedEnd
from thermoacoustic.flame import NTauFlame, CompactFlame
from thermoacoustic.solver import ThermoAcousticNetwork


# CONSTANTS

tube_length = 1
flame_position = 0.25
speed_of_sound = 343
n_flame = np.linspace(0.1, 3, 15)
tau_flame =np.linspace(0.0001, 0.015, 15)
coupling_factor = 1


inlet = ClosedEnd()
outlet = OpenEnd()
duct_one = Duct(flame_position, speed_of_sound)
duct_two = Duct(tube_length - flame_position, speed_of_sound)

omega_solutions = np.zeros((len(n_flame),len(tau_flame)), dtype=np.complex128)

i=0
j=0

for n in n_flame:
    omega_guess = 2 * np.pi * 85.75 + 1j 
    j=0
    for tau in tau_flame:
        n_tau_flame = NTauFlame(n, tau)
        compact_flame = CompactFlame(n_tau_flame, coupling_factor)

        elements = [duct_one, compact_flame, duct_two]
        network = ThermoAcousticNetwork(inlet, elements, outlet)

        omega_solution = newton(network.characteristic_value, omega_guess)
        
        omega_solutions[i][j] = omega_solution
        omega_guess = omega_solution
        j+=1
    i+=1
        
    
growth_rates = -omega_solutions.imag

N, TAU = np.meshgrid(n_flame, tau_flame, indexing='ij')

fig, ax = plt.subplots()
cf = ax.contourf(N, TAU, growth_rates, levels=80, cmap='RdBu_r')
cbar = fig.colorbar(cf, ax=ax)
cbar.set_label('Growth rate')
ax.set_xlabel('n')
ax.set_ylabel('tau [s]')
ax.set_title('Stability map (n, tau)')
ax.contour(N, TAU, growth_rates, levels=[0], colors='black', linewidths=2)
ax.text(0.15, tau_flame.max()*0.95, 'Instable (σ>0)', color='darkred', fontsize=9, va='top')
ax.text(0.15, tau_flame.max()*0.05, 'Stable (σ<0)', color='darkblue', fontsize=9, va='bottom')
fig.savefig("../results/figures/stability_map_n_tau.png", dpi=150)
plt.show()


def rayleigh_index(omega, duct_before_flame, flame):
    """Compute the sign of the Rayleigh index: Re(p* x Q') at the flame location.
    Positive = heat release in phase with pressure = destabilizing."""

    a_in, b_in = inlet.coefficients()
    if b_in != 0:
        p0, u0 = 1, -a_in / b_in
    else:
        p0, u0 = 0, 1
    state = np.array([p0, u0], dtype=complex)

    T = duct_before_flame.transfer_matrix(omega)
    state = T @ state
    p_flame, u_flame = state

    Q_prime = flame.flame_model.heat_release_response(omega) * u_flame
    RI = np.real(np.conj(p_flame) * Q_prime)
    return RI

n_test, tau_test = 1, 0.002
flame_test = CompactFlame(NTauFlame(n_test, tau_test), coupling_factor)
omega_test = newton(
    ThermoAcousticNetwork(inlet, [duct_one, flame_test, duct_two], outlet).characteristic_value,
    2 * np.pi * 85.75 + 1j
)
growth_rate_test = -omega_test.imag
RI_test = rayleigh_index(omega_test, duct_one, flame_test)

print(f"n={n_test}, tau={tau_test}: growth_rate={growth_rate_test:.2f}, Rayleigh index={RI_test:.4f}")
print(f"Signs match: {np.sign(growth_rate_test) == np.sign(RI_test)}")