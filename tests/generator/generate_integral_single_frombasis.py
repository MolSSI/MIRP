#!/usr/bin/env python3

import argparse
import random
import sys
from mpmath import mp

from common import uncontract_basis, construct_basis, all_cartesian_components, print_integral_single_input, z_map

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
parser.add_argument("--basis",    type=str, required=True, help="Path to basis set file")
parser.add_argument("--geometry", type=str, required=True, help="Path to XYZ geometry file")
parser.add_argument("--ndigits",  type=int, required=True, help="Number of digits for the value of the eri")
parser.add_argument("--ntests",   type=int, required=True, help="Number of tests to create")
parser.add_argument("--ncenter",  type=int, required=True, help="Number of centers in the integral (typically 2 or 4)")
parser.add_argument("--seed",     type=int, required=True, help="Seed to use for the pseudo-random number generator")
args = parser.parse_args()

random.seed(args.seed, version=2)

# Write out the file
with open(args.filename, 'w') as f:
    f.write("# THIS FILE IS GENERATED VIA A SCRIPT. DO NOT EDIT\n")
    f.write("#\n")
    f.write("# Input parameters for ERI generated with:\n")
    f.write("#   " + " ".join(sys.argv[:]) + "\n")
    f.write("#\n")
    f.write(str(args.ntests))
    f.write("\n")

    created_ntets = []

    with mp.workdps(args.ndigits+4):
        basis = construct_basis(args.geometry, args.basis)
        basis = uncontract_basis(basis)

        allprim = []

        for shell in basis:
            for c in all_cartesian_components(shell[1]['am']):
                allprim.append((c, shell[0], shell[1]))

        i = 0
        while i < args.ntests:
            ntet = []
            entry = []
            for n in range(args.ncenter):
                idx = random.randint(0, len(allprim)-1)
                ntet.append(idx)

                p = allprim[idx]
                Z = z_map[p[1][0].lower()]
                entry.append((Z, *p[0], *p[1][1:], p[2]['alpha'][0]))

            if not ntet in created_ntets:
                print_integral_single_input(f, entry)
                created_ntets.append(ntet)
                i += 1
