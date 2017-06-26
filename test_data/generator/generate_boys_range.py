#!/usr/bin/env python3

import argparse
import sys

import mirppy

parser = argparse.ArgumentParser()
parser.add_argument("--filename",  type=str, required=True, help="Output file name")
parser.add_argument("--max-m",     type=int, required=True, help="Maximum m value to go to")
parser.add_argument("--power",     type=int, required=True, help="Maximum power (range will be 1e-x to 1e+x)")
parser.add_argument("--ndigits",   type=int, required=True, help="Number of digits (decimal precision) required")
args = parser.parse_args()

# form the t values (as strings)
t_list = ["0"]
for i in range(-args.power, args.power+1):
    for j in range(1, 10):
        t_list.append(str(j) + "e" + str(i))

F_all = []
for t in t_list:
    F = mirppy.mirp_boys_createtest(args.max_m, t, args.ndigits)
    for m in range(0, args.max_m+1):
        F_all.append((m, t, F[m]))

if len(F_all) != ((m+1) * len(t_list)):
    raise RuntimeError("Inconsistent length of lists")

print("Calculated {} values".format(len(F_all)))

# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT")
    f.write("#\n")
    f.write("# Generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")
    f.write("#------------------------------------\n")
    f.write("# Options for {}\n".format(sys.argv[0]))
    f.write("#      Power: {}\n".format(args.power))
    f.write("#      Max m: {}\n".format(args.max_m))
    f.write("#    NDigits: {}\n".format(args.ndigits))
    f.write("#------------------------------------\n")
    f.write("#\n")
    f.write("{}\n".format(args.ndigits))

    for m, t, v in F_all:
        f.write("{} {} {}\n".format(m, t, v))
