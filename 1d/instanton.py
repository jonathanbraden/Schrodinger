import numpy as np
import matplotlib.pyplot as plt

def pot(x, lv):
    return np.cos(x) + 0.5*lv**2*np.sin(x)**2 - 1.

def phi_out(lv):
    return np.arccos((2.-lv**2)/lv**2)

if __name__=="__main__":
    pass
