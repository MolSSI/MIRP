#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--filename",  type=str, required=True, help="Output file name")
parser.add_argument("--max-m",     type=int, required=True, help="Maximum m value to go to")
parser.add_argument("--power",     type=int, required=True, help="Maximum power (range will be 1e-x to 1e+x)")
args = parser.parse_args()

# form the t values (as strings)
tlist = ["0"]
for i in range(-args.power, args.power+1):
    for j in range(1, 10):
        tlist.append(str(j) + "e" + str(i))

mtlist = []
for t in tlist:
    for m in range(0, args.max_m+1):
        mtlist.append((m, t))

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
