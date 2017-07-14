#!/usr/bin/env python3

import argparse
import sys
from mpmath import mp

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
parser.add_argument("--basis",    type=str, required=True, help="Basis set file")
parser.add_argument("--geo",      type=str, required=True, help="XYZ Geometry file")
parser.add_argument("--ndigits",  type=int, required=True, help="Number of digits for the value of the eri")
args = parser.parse_args()

mtlist = []
with mp.workdps(args.ndigits+4):
    for i in range(0, args.ntest):
        m = random.randint(0, args.max_m)
        t = random.uniform(-args.power, args.power)
        t = mp.mpf(t)    
        t = mp.power(mp.mpf(10), t)

        # min_fixed > max_fixed forces scientific notation
        t = mp.nstr(t, args.ndigits, min_fixed=1, max_fixed=0)
        mtlist.append((m,t))

# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT\n")
    f.write("#\n")
    f.write("# Values for m and t generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")

    for m,t in mtlist:
        f.write("{} {}\n".format(m, t))
