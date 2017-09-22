#!/usr/bin/env python3

import argparse
import random
import sys
from mpmath import mp

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
parser.add_argument("--max-m",    type=int, required=True, help="Maximum m value to go to")
parser.add_argument("--power",    type=int, required=True, help="Maximum power (range will be 1e-x to 1e+x)")
parser.add_argument("--seed",     type=int, required=True, help="Seed to use for the pseudo-random number generator")
parser.add_argument("--ndigits",  type=int, required=True, help="Number of digits for the value of t")
parser.add_argument("--ntests",   type=int, required=True, help="Number of tests to generate")
args = parser.parse_args()

random.seed(args.seed, version=2)

mtlist = []
with mp.workdps(args.ndigits+4):
    for i in range(0, args.ntests):
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
    f.write(str(len(mtlist)))
    f.write("\n")

    for m,t in mtlist:
        f.write("{} {}\n".format(m, t))
