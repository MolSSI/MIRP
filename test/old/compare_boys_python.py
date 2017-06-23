import mirppy
from mpmath import mp

mp.dps = 128

def boys_mpmath(n, x):
    F = []

    if x == mp.mpf("0"):
        for i in range(0, n+1):
            F.append(mp.mpf(1.0)/(mp.mpf(2.0*i+1)))
    else:
        for i in range(0, n+1):
            N = i+mp.mpf("0.5")
            F.append(mp.gammainc(N, 0, x) * 1.0/(2.0 * mp.power(x, N)))

    return F


n = 25
F = [None] * (n+1)
x = 100000.0

Fref = boys_mpmath(n, x)
mirppy.calculate_boys(F, n, x)

print("{:26.16e}    {}".format(F[n], str(Fref[n])))

