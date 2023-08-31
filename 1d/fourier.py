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
    
if __name__=="__main__":
    nLat = 1024
    lSize = 128.
    xv = np.linspace(-lSize/2.,lSize/2.,nLat, endpoint=False)

    f = np.exp(-0.5*xv**2)
    d2f = (xv**2-1.)*f
    
    mat = construct_laplacian_matrix(xv.size, 64./xv.size)
    pot = -0.5/np.cosh(xv)**2

    pot_ = BEC(xv, 1.4)
    pot = np.where(np.abs(xv)<=np.pi, pot_, 0.)
    
    op = -0.5*mat + np.diag(pot)

    w,v = np.linalg.eig(op)
    w_k, v_k = np.linalg.eig(-0.5*mat)

    ii = np.argsort(w)
    w = w[ii]
    v = v[:,ii]

    f = f / np.sqrt(np.sum(np.abs(f)**2))
    cn = v.T@f 
