import numpy as np
import matplotlib.pyplot as plt
from thermoacoustic.acoustics import Duct, ClosedEnd, OpenEnd
from thermoacoustic.solver import ThermoAcousticNetwork


# CONSTANTS
tube_length = 1
speed_of_sound = 343


inlet = ClosedEnd()
elements = [Duct(tube_length, speed_of_sound)]
outlet = OpenEnd()

network = ThermoAcousticNetwork(inlet, elements, outlet)

frequencies = np.linspace(1, 1000, 1000)

residues = []

for f in frequencies:
    omega = 2 * np.pi * f
    residues.append(np.abs(network.characteristic_value(omega)))



fig, ax = plt.subplots()
ax.plot(frequencies, residues)
ax.set_xlabel ('Frequency [Hz]')
ax.set_ylabel('|Characteristic value|')
ax.set_title('Acoustic Resonance Search')
ax.set_yscale('log')
ax.grid()

fig.savefig("../results/figures/acoustics_validation.png", dpi=150)
plt.show()


minima_indices = []

for i in range(1, len(residues) - 1):
    if residues[i] < residues[i-1] and residues[i] < residues[i+1] and residues[i] < 0.01:
        minima_indices.append(i)

found_frequencies = frequencies[minima_indices]

n_modes = len(found_frequencies)
analytical_frequencies = [(2*n - 1) * speed_of_sound / (4 * tube_length) for n in range(1, n_modes + 1)]

print("Numerical eigenfrequencies (Hz):  ", [round(f, 2) for f in found_frequencies])
print("Analytical eigenfrequencies (Hz): ", [round(f, 2) for f in analytical_frequencies])

errors = [abs(f_num - f_ana) / f_ana * 100 for f_num, f_ana in zip(found_frequencies, analytical_frequencies)]
print("Relative errors (%):              ", [round(e, 3) for e in errors])