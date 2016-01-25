import numpy as np
import cmath
import fft_real as myfft

def dct(vector):

    n = len(vector)
    fvec = np.zeros(n)

    for k in range(n):
	for i in range(n):
	    fvec[k] += vector[i]*np.cos(np.pi*k*(i+0.5)/n);

    return fvec


def fct(vector):

    n = len(vector)
    
    # auxiliary array
    aux = np.zeros(n);
    for i in range(0, n):
	aux[i] = np.sin(np.pi*(i+0.5)/n)*(vector[i] - vector[n-i-1]) + 0.5*(vector[i] + vector[n-i-1])

    vCompl = myfft.realfft(aux)

    # rearrange to give cosine transform
    for i in range(0, n/2):
	vector[2*i] = np.cos(np.pi*i/n)*vCompl[i].real + np.sin(np.pi*i/n)*vCompl[i].imag


    vector[n-1] = -0.5*vCompl[n/2].real
    i = n/2-1;
    while i > 0:
	vector[2*i-1] = vector[2*i+1] - np.sin(np.pi*i/n)*vCompl[i].real + np.cos(np.pi*i/n)*vCompl[i].imag
	vector[2*i+1] *= -1
	i -= 1

    vector[1] *= -1

    return vector


n = 2**8
omega = 2.*np.pi/n #0.001
t = np.arange(n)
x = np.cos(omega*t) #+ np.sin(omega*t) + np.sqrt(t)

discrete = dct(x)
fast = fct(x)

print len(fast), len(discrete)

print discrete[100], discrete[101], discrete[102], discrete[103]
print fast[100], fast[101], fast[102], fast[103]
print "allclose =", np.allclose(fast, discrete)

"""
for i in range(n):
    print fast[i]
"""
