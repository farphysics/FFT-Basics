import numpy as np
import time

n = 2**22
#omega = 0.001
#t = np.arange(n)
x = np.zeros(n) #np.cos(omega*t) + 1j*np.sin(omega*t)

start = time.clock()
libfft = np.fft.fft(x)
print time.clock() - start


"""
print libfft[100]

cfft = np.loadtxt("real.dat", dtype=float) + 1j*np.loadtxt("imag.dat", dtype=float)

for i in range(30):
    print libfft[i], cfft[i]

print np.allclose(libfft, cfft)

x = np.zeros(2**24)
y = np.zeros(2**24)
#x = np.random.rand(2**24)
#y = np.random.rand(2**24)
r = x+1j*y
np.fft.fft(r)
print r[100]
"""
