#!/usr/bin/env python3

import argparse
import random
import sys
from mpmath import mp
from common import print_integral_input

parser = argparse.ArgumentParser()
parser.add_argument("--filename",    type=str, required=True, help="Output file name")
parser.add_argument("--max-am",      type=int, required=True, help="Maximum AM of the basis functions")
parser.add_argument("--max-nprim",   type=int, required=True, help="Maximum number of primitives in a shell")
parser.add_argument("--max-ngen",    type=int, required=True, help="Maximum number of general contractions a shell")
parser.add_argument("--alpha-power", type=int, required=True, help="Maximum power of the exponent (range will be 1e-x to 1e+x)")
parser.add_argument("--coeff-power", type=int, required=True, help="Maximum power of the contraction coefficients (range will be 1e-x to 1e+x)")
parser.add_argument("--xyz-power",   type=int, required=True, help="Maximum power of the coordinates (range will be -1e+x to 1e+x)")
parser.add_argument("--seed",        type=int, required=True, help="Seed to use for the pseudo-random number generator")
parser.add_argument("--ndigits",     type=int, required=True, help="Number of digits for the value of the integral")
parser.add_argument("--ncenter",     type=int, required=True, help="Number of centers in the integral (typically 2 or 4)")
parser.add_argument("--ntests",      type=int, required=True, help="Number of tests to generate")
args = parser.parse_args()

random.seed(args.seed, version=2)


def generate_basis_function():
    with mp.workdps(args.ndigits+4):
        am = random.randint(0, args.max_am)

        x = mp.power(mp.mpf(10), random.uniform(-args.xyz_power, args.xyz_power))
        y = mp.power(mp.mpf(10), random.uniform(-args.xyz_power, args.xyz_power))
        z = mp.power(mp.mpf(10), random.uniform(-args.xyz_power, args.xyz_power))
        x *= mp.power(-1, random.randint(1, 2))
        y *= mp.power(-1, random.randint(1, 2))
        z *= mp.power(-1, random.randint(1, 2))

        x = mp.nstr(x, args.ndigits, min_fixed=1, max_fixed=0)
        y = mp.nstr(y, args.ndigits, min_fixed=1, max_fixed=0)
        z = mp.nstr(z, args.ndigits, min_fixed=1, max_fixed=0)

        nprim = random.randint(1, args.max_nprim)
        ngen = random.randint(1, args.max_ngen)

        alphas = []
        coeffs = []
        for i in range(nprim):
            alpha = mp.power(mp.mpf(10), random.uniform(-args.alpha_power, args.alpha_power))
            alpha = mp.nstr(alpha, args.ndigits, min_fixed=1, max_fixed=0)
            alphas.append(alpha)

            coefftmp = []
            for j in range(ngen):
                coeff = mp.power(mp.mpf(10), random.uniform(-args.coeff_power, args.coeff_power))
                coeff = mp.nstr(coeff, args.ndigits, min_fixed=1, max_fixed=0)
                coefftmp.append(coeff)
            coeffs.append(coefftmp)
                
        # Transposes the coefficient matrix
        coeffs = list(map(list, zip(*coeffs)))

        nprim = len(alphas)
        ngen = len(coeffs)


        return (am, nprim, ngen, x, y, z, alphas, coeffs)


# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT\n")
    f.write("#\n")
    f.write("# Input parameters for integral generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")
    f.write(str(args.ntests))
    f.write("\n")

    for i in range(args.ntests):
        shells = []
        for n in range(args.ncenter):
            bf = generate_basis_function()
            shells.append(bf)
        print_integral_input(f, shells)
