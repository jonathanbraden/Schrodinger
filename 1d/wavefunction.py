# Temporary script for testing
# Since this now has some useful stuff, migrate it when I'm done
# along with any data I need for the paper

import numpy as np
import matplotlib.pyplot as plt

class Wavefunction:

    def __init__(self, fName, potFile, nf=3):
        self.xv, self.pot, self.damp = np.loadtxt(potFile).T
        self.nLat = self.xv.size
        
        data = np.fromfile(fName).reshape((-1,nf,2,self.nLat))
        self.psi = data[...,0,:] + 1j*data[...,1,:]
        self.wf = self.psi[:,0]
        return

    def __str__(self):
        return "Wavefunction from 1D Schrodinger simulation"

    # Fix normalization in here
    def compute_prob_region(self, region):
        mask = ( self.xv > region[0] ) & ( self.xv < region[1] )
        return np.sum( np.abs(self.wf[:,mask])**2, axis=-1 )*self.dx

    # Fix normalization in here
    def compute_prob(self):
        return np.sum( np.abs(self.wf)**2, axis=-1)*self.dx
    
def load_wf_binary(fName, n, nf=2):
    d = np.fromfile(fName).reshape((-1,nf,2,n))
    psi = d[...,0,:] + 1j*d[...,1,:]
    return psi

def load_potential(fName):
    return np.loadtxt(fName).T

def load_log(fName):
    t, _, prob, en, en_p = np.loadtxt(fName).T
    return t, prob, en, en_p

# Functionality to code up
def compute_tunnel_point(pot):
    """
    Compute the escape point for the particle (i.e. the point on the opposite side of the barrier with equal energy to the FV)
    """
    return

def compute_prob_region(psi, xv, region):
    """
    Given the wavefunction, compute the probability that it remains trapped within some region.
    """

def density_matrix(psi):
    """
    Compute the density matrix
    """
    return

def density_matrix_diff(psi):
    """
    Compute the density matrix in mean, difference rep.
    Used by Wigner transform.
    """
    return

def wigner_function(psi):
    """
    Compute the Wigner function of psi
    """
    return
    
if __name__=="__main__":
    xv, pot, damp = load_potential('potential.dat')
    t, prob, en, en_p = load_log('log.out')
    dx = xv[1]-xv[0]
    n = xv.size
    
    nf = 3
    psi = load_wf_binary('fields.bin', n, nf)
    wf = psi[:,0]

    norm = np.sum(np.abs(wf)**2,axis=-1)*dx
                
