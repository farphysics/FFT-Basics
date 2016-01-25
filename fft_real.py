import numpy as np
import cmath
import fft_routines as myfft

def realfft(vector):

    n = len(vector)

    # split real data
    vCompl = np.zeros(n/2, dtype=complex)
    for i in range(0, n, 2):
	vCompl[i/2] = vector[i] + 1j*vector[i+1]

    # fft
    vector = myfft.bitreversalFT(vCompl)

    # add (n/2)th value = (0)th value
    tmpvec = np.zeros(n/2+1, dtype=complex)
    for i in range(0, n/2):
	tmpvec[i] = vector[i]
    tmpvec[n/2] = vector[0]

    # make vector have n values again
    vector = np.zeros(n, dtype=complex)

    # rearrange
    for i in range(0, n/2+1):
	vector[i] = (tmpvec[i] + np.conj(tmpvec[n/2-i]))/2 - 0.5j*(tmpvec[i] - np.conj(tmpvec[n/2-i]))*np.exp(-2j*np.pi*i/n)

    for i in range(1, n/2):
	vector[n-i] = np.conj(vector[i])

    return vector

n = 2**18
omega = 0.001
t = np.arange(n)
x = np.cos(omega*t)# + 1j*np.sin(omega*t)
fft = realfft(x)
