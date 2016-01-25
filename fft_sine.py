import numpy as np
import cmath
import fft_real as myfft

def dst(vector): # discrete sine transform

    n = len(vector)
    fvec = np.zeros(n);

    for k in range(1, n):
	for i in range(1, n):
	    fvec[k] += vector[i]*np.sin(np.pi*k*i/n)

    return fvec


def fst(vector): # fast sine transform

    n = len(vector)
    
    # auxiliary array
    aux = np.zeros(n);
    for i in range(1, n):
	aux[i] = np.sin(i*np.pi/n)*(vector[i] + vector[n-i]) + 0.5*(vector[i] - vector[n-i])

    vCompl = myfft.realfft(aux)

    # rearrange to give sine transform
    vector[0] = -vCompl[0].imag
    vector[1] = 0.5*vCompl[0].real
    for i in range(1, n/2):
	vector[2*i] = -vCompl[i].imag
	vector[2*i+1] = vector[2*i-1] + vCompl[i].real

    return vector


n = 2**8
omega = 2.*np.pi/n#0.001
t = np.arange(n)
x = np.sin(omega*t)# + np.cos(omega*t) + np.sqrt(t)

discrete = dst(x)
fast = fst(x)

print discrete[100], discrete[101], discrete[102], discrete[103]
print fast[100], fast[101], fast[102], fast[103]
print "allclose =", np.allclose(fast, discrete)

"""
for i in range(n):
    print fast[i]
"""
