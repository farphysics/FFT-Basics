import numpy as np
import cmath

def bitreversalFT(vector, inverse=False):
    # Returns the integer whose value is the reverse of the lowest 'bits' bits of the integer 'x'.
    def reverse(x, bits):
        y = 0
        for i in range(bits):
            y = (y << 1) | (x & 1)
            x >>= 1
        return y

    # Initialization
    n = len(vector)
    levels = 0
    while True:
        if 1 << levels == n:
            break
        elif 1 << levels > n:
            raise ValueError("Length is not a power of 2")
        else:
            levels += 1
    # Now, levels = log2(n)
    exptable = [cmath.exp((2j if inverse else -2j) * cmath.pi * i / n) for i in range(n // 2)]
    vector = [vector[reverse(i, levels)] for i in range(n)]  # Copy with bit-reversed permutation


    # Radix-2 decimation-in-time FFT
    size = 2
    while size <= n:
        halfsize = size // 2
        tablestep = n // size
        for i in range(0, n, size):
            k = 0
            for j in range(i, i + halfsize):
                temp = vector[j + halfsize] * exptable[k]
                vector[j + halfsize] = vector[j] - temp
                vector[j] += temp
		#if j == i+1:
		#    print vector[j]
                k += tablestep
        size *= 2
    return vector

"""
n = 2**18
omega = 0.001#2*np.pi/n
t = np.arange(n)
x = np.cos(omega*t)# + 1j*np.sin(omega*t)
fft = bitreversalFT(x)
#libfft = np.fft.fft(x)
#print fft[100]
"""
