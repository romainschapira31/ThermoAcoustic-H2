import cantera as ct
import numpy as np
from abc import ABC, abstractmethod
from thermoacoustic.acoustics import AcousticElement



class FlameModel(ABC):
    """Base class for a flame's heat release response model.

       Every flame model must implement heat_release_response(omega), which
      returns the (complex) factor relating the fluctuating heat release Q'
      to the upstream acoustic velocity u', i.e. Q' = heat_release_response(omega) * u'.
   """
   
    @abstractmethod
    def heat_release_response(self, omega):
        """Return the complex gain relating Q' to the upstream acoustic velocity u', at the given angular frequency omega."""
        
        pass
    
class NTauFlame(FlameModel):
    """Crocco's n-tau flame response model: the heat release fluctuation
       follows the upstream acoustic velocity with a time delay tau and a
       gain n, i.e. Q'(t) = n * u'(t - tau)."""
       
    def __init__(self, n, tau):
        """Store the interaction index n and the time delay tau."""
        
        self.n = n
        self.tau = tau
        
        
    def heat_release_response(self, omega):
        """Return n * exp(-i*omega*tau), the frequency-domain form of Crocco's law."""
        
        Q_prime = self.n * np.exp(-1j * omega * self.tau)
        
        return Q_prime
    
    
class CompactFlame(AcousticElement):
    """A compact (zero-thickness) flame element: pressure is continuous
       across the flame, while the acoustic velocity jumps in proportion to
      the fluctuating heat release, as predicted by the given flame model."""
      
    def __init__(self, flame_model, coupling_factor):
        """Store the flame's heat release response model and the coupling factor relating heat release to the velocity jump."""
        
        self.flame_model = flame_model
        self.coupling_factor = coupling_factor
        
    def transfer_matrix(self, omega):
        """Return the 2x2 transfer matrix of the compact flame: p is unchanged, u jumps by coupling_factor * heat_release_response(omega) * u_in."""
        
        T = np.array([[ 1, 0],[ 0, 1 + self.coupling_factor * self.flame_model.heat_release_response(omega)]], dtype=np.complex128)
        
        return T

class FlameCalibrator:
    """
    """
    def __init__(self, fuel, oxidizer, characteristic_distance, n=1.0,
             mechanism='gri30.yaml', T_in=298.15, P_in=101325):
        """
        """
        
        self.T_in = T_in
        self.P_in = P_in
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.n = n
        self.characteristic_distance = characteristic_distance
        self.width = 0.03
        self.refine_criteria = {"ratio": 3, "slope": 0.05, "curve": 0.05}
        self.mechanism = mechanism
        
    @classmethod
    def from_config(cls, cfg):
        """Build a calibrator from a fuel configuration (see thermoacoustic.config)."""
        fuel = cfg["fuel"]
        return cls(
            fuel=fuel["name"],
            oxidizer=fuel["oxidizer"],
            characteristic_distance=cfg["flame_model"]["characteristic_distance"],
            mechanism=fuel["mechanism"],
            T_in=fuel["T_fresh"],
            P_in=fuel["pressure"],
        )
    
    

        
    def compute_flame_velocity(self, phi):
        """
        """
        
        gas = ct.Solution(self.mechanism)
        gas.TP = self.T_in, self.P_in
        gas.set_equivalence_ratio(phi, self.fuel, self.oxidizer)
        flame = ct.FreeFlame(gas, width = self.width)
        flame.set_refine_criteria(**self.refine_criteria)
        flame.solve(loglevel=0, auto=True)
        flame_speed = flame.velocity[0]
        return flame_speed
    
    def calibrate(self, phi):
        """
        """
        flame_speed = self.compute_flame_velocity(phi)
        tau = self.characteristic_distance / flame_speed
        flame_model = NTauFlame(self.n, tau)
        return flame_model
        
        
        