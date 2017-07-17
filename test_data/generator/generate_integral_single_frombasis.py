#!/usr/bin/env python3

import argparse
import sys
from mpmath import mp

from common import construct_basis

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
parser.add_argument("--basis",    type=str, required=True, help="Path to basis set file")
parser.add_argument("--geo",      type=str, required=True, help="Path to XYZ geometry file")
parser.add_argument("--ndigits",  type=int, required=True, help="Number of digits for the value of the eri")
args = parser.parse_args()

# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT\n")
    f.write("#\n")
    f.write("# Input parameters for ERI generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")

    with mp.workdps(args.ndigits+4):
        basis = construct_basis(args.geometry, args.basis)

        for shell1 in basis:
          for shell2 in basis:
            for shell3 in basis:
              for shell4 in basis:

                alpha4 = []
                for p1 in shell1['nprim']:
                  for p2 in shell2['nprim']:
                    for p3 in shell3['nprim']:
                      for p4 in shell4['nprim']:
                        alpha4.append((p1, p2, p3, p4))

                for c1 in all_cartesian_components(shell1['am']):
                  for c2 in all_cartesian_components(shell2['am']):
                    for c3 in all_cartesian_components(shell3['am']):
                      for c4 in all_cartesian_components(shell4['am']):
                        for p in alpha4:
                            print("{} {} {} {} {} {} {}".format(c1[0], c1[1], c1[2], xyz
                                 
