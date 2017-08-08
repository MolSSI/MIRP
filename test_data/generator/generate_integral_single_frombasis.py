#!/usr/bin/env python3

import argparse
import random
import sys
from mpmath import mp

from common import uncontract_basis, construct_basis, all_cartesian_components

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
parser.add_argument("--basis",    type=str, required=True, help="Path to basis set file")
parser.add_argument("--geometry", type=str, required=True, help="Path to XYZ geometry file")
parser.add_argument("--ndigits",  type=int, required=True, help="Number of digits for the value of the eri")
parser.add_argument("--ntests",   type=int, required=True, help="Number of tests to create")
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

    created_quartets = []

    with mp.workdps(args.ndigits+4):
        basis = construct_basis(args.geometry, args.basis)
        basis = uncontract_basis(basis)

        allprim = []

        for shell in basis:
            for c in all_cartesian_components(shell[1]['am']):
                allprim.append((c, shell[0], shell[1]))


        for n in range(args.ntests):
            i = random.randint(0, len(allprim)-1)
            j = random.randint(0, len(allprim)-1)
            k = random.randint(0, len(allprim)-1)
            l = random.randint(0, len(allprim)-1)
            ijkl = (i, j, k, l)

            if ijkl in created_quartets:
                continue

            pi = allprim[i]
            pj = allprim[j]
            pk = allprim[k]
            pl = allprim[l]

            f.write("{} {} {} {} {} {} {}\n".format(*pi[0], *pi[1], pi[2]['alpha'][0]))
            f.write("{} {} {} {} {} {} {}\n".format(*pj[0], *pj[1], pj[2]['alpha'][0]))
            f.write("{} {} {} {} {} {} {}\n".format(*pk[0], *pk[1], pk[2]['alpha'][0]))
            f.write("{} {} {} {} {} {} {}\n".format(*pl[0], *pl[1], pl[2]['alpha'][0]))
            f.write("\n")

            created_quartets.append(ijkl)
