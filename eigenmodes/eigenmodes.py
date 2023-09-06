import numpy as np
import matplotlib.pyplot as plt

import fourier as spec

class QuantumEigens:

    def __init__(self, nLat, dx, pot):
        self.n = nLat
        self.dx = dx
        self.pot = pot
        
        self.lap = spec.construct_laplacian_matrix(nLat, dx)
        self.ham = -0.5*self.lap + np.diag(self.pot(xv))
        self.w, self.v = self._compute_eigens()
        
    def __str__(self):
        return "Class to extract eigenmodes"

    def _construct_hamiltonian(self):
        self.ham = -0.5*self.lap + np.diag(self.pot(xv))

    def _compute_eigens(self):
        w, v = np.linalg.eig(self.ham)
        ii = np.argsort(w)
        w = w[ii]
        v = v[:,ii]

        return w,v

    # Add symmetrization of modes
