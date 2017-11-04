import random
from mpmath import mp


def generate_basis_function(max_am, alpha_power, xyz_power, ndigits):
    with mp.workdps(ndigits+4):
        l = random.randint(0, max_am)
        m = random.randint(0, max_am)
        n = random.randint(0, max_am)

        while (l+m+n) > max_am:
            l = random.randint(0, max_am)
            m = random.randint(0, max_am)
            n = random.randint(0, max_am)

        alpha = mp.power(mp.mpf(10), random.uniform(-alpha_power, alpha_power))
        x = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        y = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        z = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        x *= mp.power(-1, random.randint(1, 2))
        y *= mp.power(-1, random.randint(1, 2))
        z *= mp.power(-1, random.randint(1, 2))

        alpha = mp.nstr(alpha, ndigits, min_fixed=1, max_fixed=0)
        x = mp.nstr(x, ndigits, min_fixed=1, max_fixed=0)
        y = mp.nstr(y, ndigits, min_fixed=1, max_fixed=0)
        z = mp.nstr(z, ndigits, min_fixed=1, max_fixed=0)

        return (l, m, n, x, y, z, alpha)


def generate_basis_shell(max_am, max_nprim, max_ngen,
                         alpha_power, coeff_power, xyz_power, ndigits):
    with mp.workdps(ndigits+4):
        am = random.randint(0, max_am)

        x = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        y = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        z = mp.power(mp.mpf(10), random.uniform(-xyz_power, xyz_power))
        x *= mp.power(-1, random.randint(1, 2))
        y *= mp.power(-1, random.randint(1, 2))
        z *= mp.power(-1, random.randint(1, 2))

        x = mp.nstr(x, ndigits, min_fixed=1, max_fixed=0)
        y = mp.nstr(y, ndigits, min_fixed=1, max_fixed=0)
        z = mp.nstr(z, ndigits, min_fixed=1, max_fixed=0)

        nprim = random.randint(1, max_nprim)
        ngen = random.randint(1, max_ngen)

        alphas = []
        coeffs = []
        for i in range(nprim):
            alpha = mp.power(mp.mpf(10), random.uniform(-alpha_power, alpha_power))
            alpha = mp.nstr(alpha, ndigits, min_fixed=1, max_fixed=0)
            alphas.append(alpha)

            coefftmp = []
            for j in range(ngen):
                coeff = mp.power(mp.mpf(10), random.uniform(-coeff_power, coeff_power))
                coeff = mp.nstr(coeff, ndigits, min_fixed=1, max_fixed=0)
                coefftmp.append(coeff)
            coeffs.append(coefftmp)
                
        # Transposes the coefficient matrix
        coeffs = list(map(list, zip(*coeffs)))

        nprim = len(alphas)
        ngen = len(coeffs)


        return (am, nprim, ngen, x, y, z, alphas, coeffs)

