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
        
         