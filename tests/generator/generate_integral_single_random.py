#!/usr/bin/env python3

import argparse
import random
import sys
from mpmath import mp
from common import print_integral_single_input
from common.randomgen import generate_basis_function

parser = argparse.ArgumentParser()
parser.add_argument("--filename",    type=str, required=True, help="Output file name")
parser.add_argument("--max-am",      type=int, required=True, help="Maximum AM of the basis functions")
parser.add_argument("--alpha-power", type=int, required=True, help="Maximum power of the exponent (range will be 1e-x to 1e+x)")
parser.add_argument("--xyz-power",   type=int, required=True, help="Maximum power of the coordinates (range will be -1e+x to 1e+x)")
parser.add_argument("--seed",        type=int, required=True, help="Seed to use for the pseudo-random number generator")
parser.add_argument("--ndigits",     type=int, required=True, help="Number of digits for the value of the integral")
parser.add_argument("--ncenter",     type=int, required=True, help="Number of centers in the integral (typically 2 or 4)")
parser.add_argument("--ntests",      type=int, required=True, help="Number of tests to generate")
args = parser.parse_args()

random.seed(args.seed, version=2)


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
        entry = []
        for n in range(args.ncenter):
            bf = generate_basis_function(args.max_am, args.alpha_power, args.xyz_power, args.ndigits)
            entry.append(bf)
        print_integral_single_input(f, entry)
