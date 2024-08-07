import numpy as np
def fft(n, m):
    ret = np.zeros([n,n], dtype='complex128')
    for i in range(n):
        for j in range(n):
            ret[i,j] = np.exp(1j*i*j*2*np.pi/n)*(2**m)
    return np.round(ret.real), np.round(ret.imag)

# print(fft(5, 6))
# print(fft(9, 5)) 

n = 15
print(n, n)
matrix = fft(n, n)
for i in range(n):
    for j in range(n):
        print(int(matrix[0][i,j]), end=' ')
    print()