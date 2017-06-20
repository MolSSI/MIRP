import ctypes
import ctypes.util
from ctypes import c_int, c_double

arb_name = ctypes.util.find_library("arbint")
arblib = ctypes.CDLL(arb_name)

d_array_3 = c_double*3

arblib.primitive_eri.argtypes = [ c_int, c_int, c_int, c_double, d_array_3,
                                  c_int, c_int, c_int, c_double, d_array_3,
                                  c_int, c_int, c_int, c_double, d_array_3,
                                  c_int, c_int, c_int, c_double, d_array_3,
                                  POINTER(c_double) ]
arblib.primitive_eri.restype = c_double

def arbint_single_eri(l1, m1, n1, alpha1, A,
                      l2, m2, n2, alpha2, B,
                      l3, m3, n3, alpha3, C,
                      l4, m4, n4, alpha4, D):

  c_A = d_array_3()
  c_B = d_array_3()
  c_C = d_array_3()
  c_D = d_array_3()

  for i in range(0, 3):
    c_A[i] = A[i]
    c_B[i] = B[i]
    c_C[i] = C[i]
    c_D[i] = D[i]

  res = float()

  return arblib.primitive_eri(l1, m1, n1, alpha1, c_A, 
                              l2, m2, n2, alpha2, c_B, 
                              l3, m3, n3, alpha3, c_C, 
                              l4, m4, n4, alpha4, c_D,
                              byref(res)) 

def arbint_boys(F, n, x):
    Ftype = c_double * (n+1)
    c_F = Ftype()

    arblib.calculate_boys.argtypes = [ Ftype, c_int, c_double ]
    arblib.calculate_boys(c_F, n, x)

    for i in range(0, n+1):
      F[i] = c_F[i]
  

x1 = 0.1009671921
y1 = 0.0713958504
z1 = 0.0
x2 = 10.0 * x1
y2 = 10.0 * y1
z2 = 10.0 * z1
x3 = 25.0 * x1
y3 = 25.0 * y1
z3 = 25.0 * z1
x4 = 150.0 * x1
y4 = 150.0 * y1
z4 = 150.0 * z1

xyz1 = ( x1, y1, z1 )
xyz2 = ( x2, y2, z2 )
xyz3 = ( x3, y3, z3 )
xyz4 = ( x4, y4, z4 )

alpha1 = 130.709311
alpha2 = 13.0709311
alpha3 = 1.30709311
alpha4 = 0.130709311


prim_eri(1, 1, 1, alpha1, xyz1, 
         1, 1, 1, alpha2, xyz2,
         1, 1, 1, alpha3, xyz3,
         1, 1, 1, alpha4, xyz4,
         a)

print("HERE: {}".format(a))

