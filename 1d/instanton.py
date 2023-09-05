import numpy as np
import matplotlib.pyplot as plt

class TunnelPotential:

    def __init__(self, x0):
        self.x0 = x0
        return

    def __str__(self):
        return "Class to compute bounce action in Quantum Mechanics"

    def pot(self, x):
        return
    
    def phi_out(self):
        return

    def phi_max(self):
        return

    def m2_fv(self, phi):
        return

    def compute_action(self):
        return
    
class BECPotential(TunnelPotential):

    def __init__(self, x0, lVal):
        super().__init__(x0)
        self.lVal = lVal

    def pot(self, x):
        v = np.cos(x/self.x0) - 1. + 0.5*self.lVal**2*np.sin(x/self.x0)**2 
        return x0**2*v
        
    def phi_out(self):
        return self.x0 * np.arccos((2.-self.lVal**2)/self.lVal**2)

    def phi_max(self):
        return

    def phi_fv(self):
        return 0._dl
    
    def m2_fv(self):
        return self.lVal**2-1.

class HertzbergPotential(TunnelPotential):

    def __init__(self, x0):
        super().__init__(x0)
        return

    def pot(self, x):
        return 0.5*x**2*(1.-0.5*(x/x0)**2) / (1.+0.5*(x/x0)**4)

    # Double check this
    def phi_out(self):
        return self.x0 * np.sqrt(np.sqrt(3.)-1.)

    def phi_max(self):
        return

    def phi_fv(self):
        return 0. 
    
    def m2_fv(self):
        return 1._dl

    
def pot(x, lv):
    return np.cos(x) + 0.5*lv**2*np.sin(x)**2 - 1.

def phi_out(lv):
    return np.arccos((2.-lv**2)/lv**2)

if __name__=="__main__":
    pass
