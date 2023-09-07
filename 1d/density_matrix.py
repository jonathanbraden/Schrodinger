# Testing code to compute density matrix and wigner
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

# TO DO:
# Peel out density matrix calculations separately from Wigner

# TO DO:
# Do more testing with periodic wavefunctions
# Clean up the one remaining loop with appropriate advanced indexing
# Double check kv sign, etc.
def wigner_from_wavefunction(psi, xv, *, periodic=False, center=True):
    # Generic setup to put in a separate subroutine
    nLat, nn = xv.size, xv.size//2 + 1 
    dx = xv[1]-xv[0]
    lSize = nLat*dx

    dk = 0.5*(2.*np.pi/lSize)
    kv = 0.5 * np.fft.fftfreq(nLat) * (2.*np.pi/dx)
    
    # This part is specific to the wavefunction input
    if periodic:
        psi_pad = np.pad(psi, (nn,nn), 'wrap')
    else:
        psi_pad = np.pad(psi, (nn,nn), 'constant')
    
    dens = np.zeros( (nLat,nn), dtype=np.complex128 )
    ind = np.arange(nn)
    
    # Remove this final loop
    for i in range(nLat):
        ii = i + nn
        dens[i,:] = psi_pad[ii-ind]*np.conj(psi_pad[ii+ind])

    wig = np.fft.irfft(dens,axis=-1) / dk  # Add size parameter for odd lattices
    
    if center:
        wig = np.fft.fftshift(wig, axes=-1)
        kv = np.fft.fftshift(kv)
    
    return wig, kv

def wigner_from_density(dens, xv, *, periodic=False, center=True):
    nLat, nn = xv.size, xv.size//2 + 1
    dx = xv[1]-xv[0]
    lSize = nLat*dx
    
    dk = 0.5*(2.*np.pi/lSize)
    kv = 0.5 * np.fft.fftfreq(nLat) * (2.*np.pi/dx)
    
    # Improve this with appropriate indexing
    if periodic:
        method = 'wrap'
    else:
        method='constant'                      
    dens_pad = np.pad(dens, ((nn,nn),(nn,nn)), method)
    ind = np.arange(nn)
    dens_ = np.zeros( (nLat,nn), dtype=np.complex128 )
    for i in range(nLat):
        ii = i+nn
        dens_[i,:] = dens_pad[ii-ind,ii+ind]
    
    wig = np.fft.irfft(dens_, axis=-1) / dk
    
    if center:
        wig = np.fft.fftshift(wig, axes=-1)
        kv = np.fft.fftshift(kv)
    
    return wig, kv

if __name__=="__main__":
    n = 4
    xv = np.linspace(-0.5,0.5,n,endpoint=False)
    psi = np.sin(2.*np.pi*xv)
    
     
