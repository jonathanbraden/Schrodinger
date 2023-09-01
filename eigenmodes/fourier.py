import numpy as np
import matplotlib.pyplot as plt

def fourier_derivative_real(f, dx):

    n = f.size; norm = 2.*np.pi/dx
    fk = 1j*np.fft.rfft(f)*np.fft.rfftfreq(n)*norm
    return np.fft.irfft(fk)

def fourier_laplacian_real(f, dx):

    n = f.size; norm = 2.*np.pi/dx
    fk = -np.fft.rfft(f)*(np.fft.rfftfreq(n)*norm)**2
    return np.fft.irfft(fk)

def construct_laplacian_matrix(nLat, dx):

    norm = 2.*np.pi / dx
    kv = np.fft.rfftfreq(nLat)
    Fk = -kv**2*np.fft.rfft( np.eye(nLat), axis=-1 )
    return norm**2*np.fft.irfft(Fk, n=nLat, axis=-1)

def posch_teller(x, lv):
    return -0.5*lv*(lv+1)/np.cosh(x)**2

def Hertzberg(x, x0):
    return 0.5*x**2*(1.-0.5*(x/x0)**2)/(1.+0.5*(x/x0)**4) + 0.5*x0**2

def BEC(x, lv):
    return np.cos(x) + 1. + 0.5*lv**2*np.sin(x)**2

def eigenvalue_stuff(op):
    w, v = np.linalg.eig(op)
    
    ii = np.argsort(w)
    w = w[ii]
    v = v[:,ii]

    f = np.zeros_like(w)  # Make this the initial wavefunction

    # Figure out how to separate even an odd modes

    # Note: the c_n signs bounce around, but seems to be related to the sign
    #  of the eigenmodes at the origin.  Simply arrange them to all have the same sign (and figure out what to do with the odd modes

# Debug the ordering of the conjugate here
def density_matrix(psi):
    return np.outer(psi,np.conjugate(psi))

# One option to get the Wigner is to compute the density matrix
# as above, then do a bunch of shifting (via matrix multiplication or
# other) to put it into the mean diff form
#
# Need to check FFT conventions to make sure everything is correct

from scipy import linalg
# Fix this
# Add option for localized instead of periodic version
def density_matrix_wigner(psi, *, periodic=True):
    """
    Computes density matrix in mean field different form suitable for obtaining Wigner function
    """
    n = psi.size

    # Ugly looping version that is slow in python
    dens = np.zeros((n,n))
    # Periodic version
    #ind = np.concatenate([np.arange(n)]*3)

    #psi = np.pad(psi,(n,n), 'constant')
    psi = np.pad(psi,(n,n), 'wrap')
    psi[:n] = 0.
    psi[2*n:] = 0.
    
    for i in range(n):
        for j in range(n):
            dens[i,j] = psi[i+n-j]*np.conjugate(psi[i+n+j]) 
            
    #bins = np.arange(n)
    #indices = linalg.hankel(bins, bins+n-(n%2))

    #pad_psi = np.pad(psi, (n,n), 'constant')  # Is constant correct?
    #density_matrix = pad_psi[indices+n]*np.conjugate(pad_psi[indices[::,::-1]])
    
    return dens

# Fix this (taken from some online repo that looks buggy as hell)
def compute_wigner(psi):
    n = psi.size
    bins = np.arange(n)
    indices = linalg.hankel( bins, bins+n-(n%2) )

    pad_psi = np.pad(psi,(n,n),'constant')
    dens = pad_psi[indices+n]*np.conjugate(pad_psi[indices[::,::-1]])

    wigner = np.real(np.fft.fft(dens, axis=1)).T
    # Should probably flip frequencies
    return wigner
    
if __name__=="__main__":
    nLat = 512
    lSize = 16.*np.pi
    xv = np.linspace(-lSize/2.,lSize/2.,nLat, endpoint=False)

    # Better to normalize these differently?
    f = np.exp(-0.5*xv**2/np.sqrt(1.4**2-1.))
    f = f / np.sqrt(np.sum(np.abs(f)**2))

    # Construct Hamiltonian
    mat = construct_laplacian_matrix(xv.size, lSize/xv.size)
    pot = BEC(xv, 1.4)
    pot = np.where(np.abs(xv)<=np.pi, pot, 0.)
    op = -0.5*mat + np.diag(pot)

    # Combine these into a function
    w_free, v_free = np.linalg.eig(-0.5*mat)
    ii = np.argsort(w_free)
    w_free = w_free[ii]
    v_free = v_free[:,ii]

    w,v = np.linalg.eig(op)
    ii = np.argsort(w)
    w = w[ii]
    v = v[:,ii]

    # Do the symmetrization, etc. here
    v = v*np.sign( v[nLat//2,:] )
    f = f / np.sqrt(np.sum(np.abs(f)**2))
    cn = v.T@f

    # Check if this is the best way to do this
    mask = (np.abs(cn) > 1.e-3)

    for s2 in [0.5,0.75,1.,1.5,2.]:
        f_ = np.exp(-0.5*xv**2/s2)
        f_ = f_ / np.sqrt(np.sum(np.abs(f_)**2))
        cn_ = v.T@f_
        plt.plot(cn[mask],cn_[mask])
