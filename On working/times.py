import ctypes as ct
import numpy as np

N = 2

fortlib = ct.cdll.LoadLibrary("./flib.so")
f = fortlib.mult

f.argtypes=[ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.c_int]
f.restype = None


a = np.array(np.random.random((N,N)), order="F")
a_ptr = a.ctypes.data_as(ct.POINTER(ct.c_double))

b = np.array(np.random.random((N,N)), order="F")
b_ptr = b.ctypes.data_as(ct.POINTER(ct.c_double))

c = np.zeros((N,N), order="F")
c_ptr = c.ctypes.data_as(ct.POINTER(ct.c_double))

f(a_ptr,b_ptr,c_ptr, ct.c_int(N))

print("La matriz a:")
print(a)
print("La matriz b:")
print(b)
print("El resultado de operar a*b:")
print(c)