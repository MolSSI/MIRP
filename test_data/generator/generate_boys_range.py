#!/usr/bin/env python3

import argparse
import sys

from common import boys
from mpmath import mp

# Note that types are being stored as a string. This is so mpmath
# can parse it without turning it into a (possibly) lesser-precision float
parser = argparse.ArgumentParser()
parser.add_argument("--filename",     type=str, required=True, help="Output file name")
parser.add_argument("--max-m",        type=int, required=True, help="Maximum m value to go to")
parser.add_argument("--power",        type=int, required=True, help="Maximum power (range will be 1e-x to 1e+x)")
parser.add_argument("--prec",         type=int, required=True, help="Final precision to truncate to (in bits)")
parser.add_argument("--working-prec", type=int, required=True, help="Binary precision (in bits)")
args = parser.parse_args()

# Set the precision of the math library
mp.prec = args.prec

# List of m values
m_list = list(range(0, args.max_m+1))

# form the t values (as strings)
t_list = [mp.mpf("0")]
for i in range(-args.power, args.power+1):
    for j in range(1, 10):
        t_list.append(mp.mpf(str(j) + "e" + str(i)))

# Temporarily work in working precision
with mp.workprec(args.working_prec):
    F_list = boys.calculate_boys_list(m_list, t_list)
    if len(F_list) != (len(m_list) * len(t_list)):
        raise RuntimeError("Inconsistent length of list returned from calculate_boys_list")

print("Calculated {} values".format(len(F_list)))

# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT")
    f.write("#\n")
    f.write("# Generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")
    f.write("#------------------------------------\n")
    f.write("# Options for {}\n".format(sys.argv[0]))
    f.write("#              Power: {}\n".format(args.power))
    f.write("#              Max m: {}\n".format(args.max_m))
    f.write("#          Precision: {}\n".format(args.prec))
    f.write("#  Working Precision: {} bits\n".format(args.working_prec))
    f.write("#------------------------------------\n")
    f.write("#\n")
    f.write("{}\n".format(args.prec))

    for m, t, v in F_list:
        f.write("{} {} {}\n".format(m,
                                    mp.nstr(t, mp.dps, max_fixed = 1, min_fixed = 2, strip_zeros=False),
                                    mp.nstr(v, mp.dps, max_fixed = 1, min_fixed = 2, strip_zeros=False)))
