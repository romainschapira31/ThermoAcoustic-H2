import cantera as ct
import numpy as np

gas = ct.Solution('gri30.yaml')

T = 300
P = ct.one_atm

gas.set_equivalence_ratio(phi=1.0, fuel='CH4', oxidizer='O2:1.0, N2:3.76')

gas.TP = T, P



flame = ct.FreeFlame(gas, width=0.03)

flame.solve(loglevel=1, auto=True)

speed_flame = flame.velocity[0]
max_T = max(flame.T)

print(f"Vitesse de la flamme : {speed_flame * 100:.2f} cm/s")
print(f"Température des gaz brûlés : {max_T:.0f} K")