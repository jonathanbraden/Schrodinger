# Temporary script for testing
# Since this now has some useful stuff, migrate it when I'm done
# along with any data I need for the paper

import numpy as np
import matplotlib.pyplot as plt

# Debug these
def laplacian_spectral(f, dx, *, axis=1):
    n = f.shape[axis]
    
    norm = (2.*np.pi / dx)
    kv = np.fft.fftfreq(n) * norm
    return np.fft.ifft( -kv**2 * np.fft.fft(f, axis=axis), axis=axis)

def grad_spectral(f, dx, *, axis=1):
    n = f.shape[axis]
    norm = (2.*np.pi / dx)
    kv = np.fft.fftfreq(n) * norm
    return np.fft.ifft( 1j*kv * np.fft.fft(f, axis=axis), axis=axis)
    return

# Add real gradient for current divergence



class Wavefunction:

    def __init__(self, fName, potFile, logFile, nf=3):
        self.xv, self.pot, self.damp = np.loadtxt(potFile).T
        self.nLat = self.xv.size
        self.dx = self.xv[1]-self.xv[0]

        # Simplify this at some point
        self.t, self.prob, self.en, self.en_p = load_log('log_asym.out')

        data = np.fromfile(fName).reshape((-1,nf,2,self.nLat))
        self.psi = data[...,0,:] + 1j*data[...,1,:]
        self.wf = self.psi[:,0]
        return

    def __str__(self):
        return "Wavefunction from 1D Schrodinger simulation"

    def compute_prob_region(self, region):
        mask = ( self.xv > region[0] ) & ( self.xv < region[1] )
        return np.sum( np.abs(self.wf[:,mask])**2, axis=-1 )*self.dx

    def compute_prob(self):
        return np.sum( np.abs(self.wf)**2, axis=-1)*self.dx

    def momentum_operator(self):
        """
        p = -i\hbar\partial_x

        which is Im(\partial\psi)
        """
        return -1j*derivative_spectral(self.wf, self.dx)
    
    # Check ordering on this
    def density_matrix(self):
        return np.outer(self.wf, np.conj(self.wf))

    # Debug normalization, and check sign convention
    def prob_current(self):
        """
        Returns the probability current
         J = -\hbar/(2mi)(\psi^*\nabla\psi - \psi\nabla\psi^*) = 
        """
        grad = grad_spectral(self.wf, self.dx)
        return np.imag(np.conj(self.wf)*grad)

    # Debug normalization and check sign convention
    def prob_current_divergence(self, *, method='laplacian'):
        if method=='laplacian':
            lap = laplacian_spectral(self.wf, self.dx)
            div = -np.imag(np.conj(self.wf)*lap)
        # Fix this to use the real gradient
        elif method=='divergence':
            grad = self.prob_current()
            div = -grad_spectral(grad, self.dx) # Oops, real derivative here
        else:
            print("Invalid method")
            return None
        return div

    def overlap_ic(self):
        """
        Compute overlap with initial state
        """
        return (self.wf @ self.wf[0])*dx

    def overlap_w_psi(self, psi):
        return (self.wf @ psi)*dx
    
    # Add a bunch of plotting routines


    
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
    wf_sym = Wavefunction('fv_sym.bin', 'pot_sym.dat', 'log_sym.out')
    wf_asym = Wavefunction('fv_asym.bin', 'pot_asym.dat', 'log_asym.out')       
